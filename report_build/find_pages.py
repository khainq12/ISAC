import fitz
import io, sys
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
doc = fitz.open('LyThuyet_AFDM_OTFS_OFDM.pdf')
print('total pages', doc.page_count)
for i, page in enumerate(doc):
    text = page.get_text()
    if 'Cam bien' in text or 'CRB' in text or 'Muc 13' in text or 'Cảm biến' in text:
        print(i+1, text[:100].replace(chr(10),' | '))
