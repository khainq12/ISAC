import fitz
import sys
doc = fitz.open('LyThuyet_AFDM_OTFS_OFDM.pdf')
pages = [int(x) for x in sys.argv[1:]]
for p in pages:
    page = doc[p-1]
    pix = page.get_pixmap(dpi=150)
    pix.save(f'page-{p:02d}.png')
print('done')
