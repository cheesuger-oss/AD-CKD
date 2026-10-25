"""Reuse the original Python generator; disable only colourbar rasterization."""
from pathlib import Path
import runpy
from matplotlib.colorbar import Colorbar

audit = Path(__file__).resolve().parent
original = Colorbar._draw_all

def vector_draw(self):
    original(self)
    if self.solids is not None:
        self.solids.set_rasterized(False)
        # Avoid PDF-viewer hairline seams between adjacent vector colour cells.
        self.solids.set_edgecolor('face')

Colorbar._draw_all = vector_draw
ns = runpy.run_path(str(audit / 'assemble_supplementary_figure_s3_AAD.py'))
stage = audit / 'vector_export_stage'
stage.mkdir(exist_ok=True)
ns['main'].__globals__.update(ROOT=stage,
    OUT_PDF=stage/'Supplementary_Figure_S3_animal_enrichment_summary.pdf',
    OUT_PNG=stage/'Supplementary_Figure_S3_animal_enrichment_summary.png')
ns['main']()
