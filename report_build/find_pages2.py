import fitz, io, sys
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
doc = fitz.open('LyThuyet_AFDM_OTFS_OFDM.pdf')
print('total pages', doc.page_count)
markers = ['Vì sao cần', 'Định nghĩa bậc phân tập', 'AFDM đạt bậc phân tập', 'OTFS: bậc phân tập giảm',
           'Cấu trúc ma trận kênh hiệu dụng', 'BER/SER trên một kênh đại diện', 'Bậc phân tập qua thu ML',
           'Cảm biến: sai số ước lượng', 'Ảnh hưởng của tốc độ mục tiêu', 'Giới thiệu mô phỏng', 'Kết luận']
for i, page in enumerate(doc):
    text = page.get_text()
    for m in markers:
        if m in text:
            print(i+1, m)
