import sys
import os

def _append_relative_path(path):
    sys.path.append(os.path.join(
        os.path.dirname(os.path.realpath(__file__)), path))

try:
    import wget_provider
except ImportError:
    _append_relative_path('wget_provider')
    import wget_provider

try:
    import copernicus_download
except ImportError:
    _append_relative_path('copernicus_download')
    import copernicus_download
