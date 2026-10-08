import fitz, io, sys
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
doc = fitz.open('LyThuyet_AFDM_OTFS_OFDM.pdf')
print('total pages', doc.page_count)
markers = ['AFDM đạt bậc phân tập đầy đủ', 'OTFS: bậc phân tập giảm', 'Tham số và phương pháp mô phỏng',
           'Sơ đồ tổng quan', 'Ảnh hưởng của tốc độ mục tiêu', 'Thảo luận', 'Kết luận',
           'Phụ lục A']
for i, page in enumerate(doc):
    text = page.get_text()
    for m in markers:
        if m in text:
            print(i+1, m)
