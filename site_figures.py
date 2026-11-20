"""site_figures.py -- the website's figure files: WebP, 800 pixels wide (the site must stay under the 1 GB Pages limit)."""
import os

MAX_WIDTH = 800        # pixels
QUALITY = 75           # WebP quality: 54 kB per transit figure, 30 kB per EB figure (JPEG 80: 103 and 68 kB)
TRUE_TIME = 'true_time_fits.tsv'   # rows listed at a fit on the samples' true times (true_time_catalog.py table)
_true_time = {}


def to_webp(source, target, max_width=MAX_WIDTH, quality=QUALITY):
    """Write the image `source` as WebP at `target`, shrunk to `max_width` pixels."""
    from PIL import Image
    os.makedirs(os.path.dirname(target), exist_ok=True)
    image = Image.open(source).convert('RGB')
    if image.width > max_width:
        image = image.resize((max_width, int(image.height * max_width / image.width)), Image.LANCZOS)
    image.save(target + '.tmp', 'WEBP', quality=quality)
    os.replace(target + '.tmp', target)


def true_time_plots(here):
    """The figures (WebP names) of the rows listed at a true-time fit."""
    if here not in _true_time:
        path = os.path.join(here, TRUE_TIME)
        plots = set()
        if os.path.exists(path):
            import pandas as pd
            rows = pd.read_csv(path, sep='\t')
            plots = set(rows.loc[rows['adopted'], 'plot'].map(webp))
        _true_time[here] = plots
    return _true_time[here]


def install(source, here, plot, true_time=False):
    """Publish the image `source` as the WebP of the row figure `plot` (a path under the site directory `here`) and
    remove the file of the same name in another format that it replaces. The figure of a row listed at a true-time
    fit is replaced only by a fold on the true times (`true_time`): a grid fold is offset from that ephemeris.
    Returns whether the figure was written."""
    if not true_time and webp(plot) in true_time_plots(here):
        return False
    to_webp(source, os.path.join(here, webp(plot)))
    old = os.path.join(here, plot)
    if old != os.path.join(here, webp(plot)) and os.path.exists(old):
        os.remove(old)
    return True


def webp(plot):
    """The WebP name of a figure path."""
    return os.path.splitext(plot)[0] + '.webp'


def published(here, plot):
    """The figure of a row as the site should link it: the WebP once it exists, else the file named, else ''."""
    if not isinstance(plot, str) or plot == '':
        return ''
    for name in (webp(plot), plot):
        if os.path.exists(os.path.join(here, name)):
            return name
    return ''
