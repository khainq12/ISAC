import fitz, io, sys
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
doc = fitz.open('LyThuyet_AFDM_OTFS_OFDM.pdf')
markers = ['Giới thiệu mô phỏng', 'Tham số mô phỏng tổng hợp', 'Phương pháp đối chứng',
           'Cấu trúc ma trận kênh hiệu dụng', 'BER/SER trên một kênh đại diện',
           'Cảm biến: sai số ước lượng', 'Ảnh hưởng của tốc độ mục tiêu',
           'Thảo luận', 'Hạn chế của nghiên cứu', 'Kết luận']
for i, page in enumerate(doc):
    text = page.get_text()
    for m in markers:
        if m in text:
            print(i+1, m)
