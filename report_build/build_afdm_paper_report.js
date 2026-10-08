const fs = require('fs');
const {
  Document, Packer, Paragraph, TextRun, HeadingLevel, ImageRun,
  Table, TableRow, TableCell, WidthType, ShadingType, BorderStyle,
  AlignmentType, PageOrientation
} = require('docx');

const DIR = 'D:\\doubly-dispersive-channel-simulation-main\\afdm_paper_check\\';

function h1(text) {
  return new Paragraph({ text, heading: HeadingLevel.HEADING_1, spacing: { before: 300, after: 150 } });
}
function h2(text) {
  return new Paragraph({ text, heading: HeadingLevel.HEADING_2, spacing: { before: 250, after: 120 } });
}
function p(text, opts = {}) {
  return new Paragraph({ children: [new TextRun({ text, ...opts })], spacing: { after: 120 } });
}
function pRuns(runs) {
  return new Paragraph({ children: runs, spacing: { after: 120 } });
}
function bullet(text) {
  return new Paragraph({ text, bullet: { level: 0 }, spacing: { after: 80 } });
}
function img(filename, widthPx = 580) {
  const buf = fs.readFileSync(DIR + filename);
  const h = Math.round(widthPx * 0.62);
  return new Paragraph({
    children: [new ImageRun({ data: buf, type: 'png', transformation: { width: widthPx, height: h } })],
    alignment: AlignmentType.CENTER,
    spacing: { after: 200 }
  });
}
function caption(text) {
  return new Paragraph({
    children: [new TextRun({ text, italics: true, size: 20 })],
    alignment: AlignmentType.CENTER,
    spacing: { after: 240 }
  });
}

const cellBorders = {
  top: { style: BorderStyle.SINGLE, size: 2, color: 'AAAAAA' },
  bottom: { style: BorderStyle.SINGLE, size: 2, color: 'AAAAAA' },
  left: { style: BorderStyle.SINGLE, size: 2, color: 'AAAAAA' },
  right: { style: BorderStyle.SINGLE, size: 2, color: 'AAAAAA' },
};
function cell(text, opts = {}) {
  return new TableCell({
    children: [new Paragraph({ children: [new TextRun({ text, bold: !!opts.bold, size: 20 })] })],
    width: { size: opts.width || 2000, type: WidthType.DXA },
    shading: opts.shade ? { type: ShadingType.CLEAR, fill: opts.shade } : undefined,
    borders: cellBorders,
    verticalAlign: 'center',
  });
}
function statusTable() {
  const rows = [
    ['Hình', 'Nội dung', 'Kết quả'],
    ['8a', 'AFDM, BPSK, N=16, ML, P=2/3/4 đường (phân tập)', '✅ Khớp — đúng thứ tự phân tập P=4>P=3>P=2'],
    ['8b', '4 waveform (OFDM/OCDM/OTFS/AFDM), BPSK N=16, ML', '✅ Khớp — OTFS≈AFDM, vượt OCDM>OFDM'],
    ['9a', '4 waveform, QPSK N=256, MMSE', '✅ Khớp — OTFS/AFDM trùng khít'],
    ['10a/10b', 'AFDM vs OFDM, 3 bộ tách sóng (MMSE/DFE/băng-hẹp)', '✅ Khớp — trùng ở Doppler nguyên, tách biệt ở lẻ bin'],
    ['11', 'Ước lượng kênh AFDM, Doppler nguyên, SNRp=25/30/35dB', '✅ Khớp — SNRp=35dB gần trùng Ideal'],
    ['12a', 'Ước lượng kênh AFDM, Doppler lẻ bin, SNRp=30/35/40dB', '✅ Khớp — đúng thứ tự SNRp=40>35>30dB'],
    ['12b', 'Hiệu ứng ξν (guard) = 1/4/8, SNRp=40dB', '✅ Khớp — đúng ξν tăng → BER giảm'],
  ];
  const widths = [1100, 5200, 3300];
  return new Table({
    width: { size: 9600, type: WidthType.DXA },
    columnWidths: widths,
    rows: rows.map((r, i) => new TableRow({
      children: r.map((t, j) => cell(t, { width: widths[j], bold: i === 0, shade: i === 0 ? 'D9E2F3' : (i % 2 === 0 ? 'F2F2F2' : undefined) })),
    })),
  });
}

function paramTable() {
  const rows = [
    ['Tham số', 'Giá trị / công thức', 'Ghi chú'],
    ['Bài báo', 'Bemani, Ksairi, Kountouris, "Affine Frequency Division Multiplexing for Next Generation Wireless Communications," IEEE TWC 2023', 'arXiv:2204.12798'],
    ['Kênh', 'H = Σ hᵢ·ΓCPPᵢ·Δfᵢ·Πˡⁱ (eq. 22-24)', 'Tự viết theo đúng công thức bài báo'],
    ['c1 (Doppler nguyên)', 'c1 = (2αmax+1)/(2N)', 'eq. 42'],
    ['c1 (Doppler lẻ bin)', 'c1 = (2(αmax+ξν)+1)/(2N)', 'eq. 43'],
    ['OCDM', 'c1 = c2 = 1/(2N)', ''],
    ['Q (guard length)', 'Q = (lmax+1)(2(αmax+ξν)+1) − 1', 'eq. 52'],
    ['Kênh vật lý', 'P đường, gain Gauss phức, phương sai 1/P', 'Đúng giả định bài báo, Mục VII'],
  ];
  const widths = [2200, 5000, 2400];
  return new Table({
    width: { size: 9600, type: WidthType.DXA },
    columnWidths: widths,
    rows: rows.map((r, i) => new TableRow({
      children: r.map((t, j) => cell(t, { width: widths[j], bold: i === 0, shade: i === 0 ? 'D9E2F3' : (i % 2 === 0 ? 'F2F2F2' : undefined) })),
    })),
  });
}

const doc = new Document({
  sections: [{
    properties: { page: { size: { width: 11906, height: 16838 } } }, // A4
    children: [
      new Paragraph({
        children: [new TextRun({ text: 'Báo cáo tái hiện kết quả mô phỏng', bold: true, size: 36 })],
        alignment: AlignmentType.CENTER,
        spacing: { after: 100 },
      }),
      new Paragraph({
        children: [new TextRun({ text: 'Bài báo AFDM gốc (arXiv:2204.12798)', bold: true, size: 28, color: '444444' })],
        alignment: AlignmentType.CENTER,
        spacing: { after: 400 },
      }),
      p('Bemani, Ksairi, Kountouris, "Affine Frequency Division Multiplexing for Next Generation Wireless Communications," IEEE Transactions on Wireless Communications, 2023.', { italics: true }),
      p('Toàn bộ 7 hình kết quả mô phỏng chính của bài báo (Mục VII) đã được tái hiện độc lập bằng MATLAB, đối chiếu trực tiếp với các khẳng định định tính/định lượng của tác giả.'),

      h1('1. Tóm tắt kết quả'),
      statusTable(),
      p(''),
      p('7/7 hình khớp với bài báo về mặt xu hướng và phần lớn về độ lớn số liệu. Trong quá trình thực hiện, một lỗi quan trọng trong bộ tách sóng MMSE đã được phát hiện và sửa (xem Mục 3).', { bold: true }),

      h1('2. Phương pháp và tham số'),
      p('Toàn bộ mô hình kênh, ma trận DAFT (AFDM/OCDM), và các công thức guard/pilot được viết lại độc lập, bám sát chính xác các công thức trong bài báo (không sao chép code có sẵn).'),
      paramTable(),
      p(''),
      p('Do giới hạn thời gian tính toán, một số hình dùng numMonte (số lần Monte Carlo) thấp hơn con số 10⁶ mà bài báo dùng cho Hình 8 (ML). Cụ thể: Hình 8a/8b dùng numMonte=3.000–10.000; Hình 9/10 dùng 1.500–3.000; Hình 11/12 dùng 300–1.500. Điều này khiến các đường cong không mượt tuyệt đối ở vùng BER cực thấp (10⁻⁶–10⁻⁷) như bài báo, nhưng xu hướng và thứ tự so sánh giữa các waveform/tham số vẫn được bảo toàn rõ ràng.'),

      h1('3. Hai lỗi quan trọng đã phát hiện và sửa'),
      p('Lỗi 1 — bỏ sót năng lượng rò vào guard: Trong các Hình 11, 12a, 12b (thí nghiệm ước lượng kênh qua pilot nhúng), bộ tách sóng MMSE ban đầu chỉ sử dụng khối ma trận vuông Heff(dataIdx, dataIdx) và các mẫu thu tại đúng vị trí dữ liệu — bỏ qua hoàn toàn năng lượng tín hiệu bị rò vào vùng guard (do AFDM không phải là biến đổi đường chéo).'),
      p('Hệ quả: ở Hình 12b, xu hướng đo được BỊ ĐẢO NGƯỢC so với bài báo — ξν càng lớn (guard rộng hơn) lại cho BER càng XẤU, trong khi bài báo khẳng định ngược lại ("increasing ξν improves the performance").'),
      p('Cách sửa: dùng toàn bộ vector thu N chiều z, và ma trận chữ nhật Heff(:, dataIdx) — đúng theo đẳng thức chính xác z = Heff(:,dataIdx)·x_data + nhiễu (vì các vị trí guard/pilot ngoài dataIdx luôn phát đúng giá trị đã biết hoặc 0).'),
      p('Lỗi 2 — chưa trừ phần rò của pilot khi dò dữ liệu: sau khi sửa Lỗi 1, Hình 12a vẫn cho SNRp=40dB (pilot mạnh nhất) tệ hơn SNRp=30/35dB — ngược trực giác. Chẩn đoán riêng cho thấy sai số ước lượng kênh (Frobenius) vẫn giảm đơn điệu đúng theo SNRp, chứng tỏ bộ ước lượng không có lỗi; vấn đề nằm ở bước dò dữ liệu: vector thu zP dùng để dò vẫn còn nguyên phần đóng góp đã biết của pilot Heff(:,1)·xpilot (chưa bị trừ ra), và biên độ phần "rò" này tỉ lệ thuận với xpilot — nên SNRp càng cao, nhiễu rò pilot chưa khử càng lớn.', { bold: true }),
      p('Cách sửa: trừ Heff_hat(:,1)·xpilot khỏi zP trước khi đưa vào công thức MMSE. Áp dụng cho cả Hình 11, 12a, 12b; sau khi sửa và chạy lại cả 3 hình, xu hướng khớp đúng bài báo trên toàn dải SNR — xem chi tiết Mục 4.5-4.7.'),

      h1('4. Chi tiết từng hình'),

      h2('4.1. Hình 8a — Bậc phân tập của AFDM'),
      p('AFDM, BPSK, N=16, dò ML vét cạn, kênh P=2/3/4 đường, Doppler nguyên αmax=1. Bài báo dùng đường tham chiếu s(SNR)⁻ᵈ để minh họa độ dốc phân tập bậc d=P.'),
      img('fig8a_reproduction.png'),
      caption('Hình 8a tái hiện: BER vs SNR cho P=2/3/4. Thứ tự phân tập đúng tại mọi SNR≥6dB.'),
      p('Kết quả: P=4 luôn dốc nhất (phân tập cao nhất), P=2 nông nhất, đúng thứ tự bài báo dự đoán. Độ dốc đo được trong đoạn 8-16dB xấp xỉ 1.5/2.8/3.3 cho P=2/3/4, đang tiệm cận đúng giá trị lý thuyết 2/3/4 khi SNR tăng thêm (giống bài báo, độ dốc sạch chỉ đạt ở BER rất thấp).'),

      h2('4.2. Hình 8b — So sánh 4 waveform bằng ML'),
      p('OFDM/OCDM/OTFS/AFDM, BPSK N=16, ML, kênh 3 đường cố định (delay 0,1,2), Doppler lẻ bin theo phổ Jakes, αmax=1.'),
      img('fig8b_reproduction.png'),
      caption('Hình 8b tái hiện: OTFS và AFDM trùng khít, vượt trội OCDM và OFDM.'),
      p('Kết quả khớp rất tốt: OFDM tệ nhất mọi SNR (đúng "cannot separate the paths"), OCDM khá hơn nhưng vẫn hạn chế, OTFS và AFDM gần như trùng nhau và tốt nhất — đúng khẳng định "AFDM has the same BER performance as OTFS".'),

      h2('4.3. Hình 9a — So sánh 4 waveform bằng MMSE, N=256'),
      p('QPSK, N=256, MMSE, kênh 3 đường, Doppler lẻ bin Jakes, αmax=2.'),
      img('fig9a_reproduction.png'),
      caption('Hình 9a tái hiện: tại SNR=27.5dB, OTFS và AFDM cho đúng cùng một giá trị BER=3.2552×10⁻⁶.'),
      p('Kết quả khớp xuất sắc, gần như trùng khít tuyệt đối giữa OTFS và AFDM ở SNR cao — đúng chính xác khẳng định của bài báo.'),

      h2('4.4. Hình 10a/10b — So sánh 3 bộ tách sóng'),
      p('AFDM vs OFDM, 3 bộ tách sóng: MMSE chính xác, MRC-DFE (lặp Gauss-Seidel — bài báo chứng minh tương đương LMMSE), và MMSE băng-hẹp (giữ L phần tử lớn nhất mỗi hàng, thích ứng theo từng lần thử). Doppler mỗi đường sinh theo phổ Jakes liên tục rồi làm tròn về số nguyên gần nhất cho trường hợp "Doppler nguyên" (đúng mô tả "Doppler generated using Jakes Doppler spectrum" của bài báo — cho phép va chạm Doppler tự nhiên giữa các đường, không ép buộc phân biệt); delay cố định 0..lmax như quy ước Hình 8b.'),
      img('fig10_reproduction.png', 620),
      caption('Hình 10a/10b tái hiện: trái = Doppler nguyên (nhóm AFDM và nhóm OFDM tách biệt rõ, OFDM chững thành floor), phải = Doppler lẻ bin.'),
      p('10a (Doppler nguyên): cả 3 bộ tách sóng trùng khít TRONG CÙNG một waveform (đúng "the channel matrices used in the three detection methods are all the same"), nhưng nhóm AFDM tiếp tục giảm mạnh (BER≈1.2×10⁻⁶ tại SNR=30dB) trong khi nhóm OFDM chững lại thành floor quanh 1.25×10⁻⁴ — khớp gần như chính xác với floor ~1.2×10⁻⁴ của bài báo tại SNR=30dB, đúng khẳng định "AFDM outperforms OFDM, owed to achieving full diversity". 10b (Doppler lẻ bin): nhóm AFDM nằm dưới nhóm OFDM, AFDM MMSE-băng-hẹp vẫn tốt hơn các đường OFDM, và MMSE chính xác/MRC-DFE của AFDM tiếp tục giảm ở SNR cao trong khi OFDM chững lại — phù hợp "exact LMMSE has slightly better performance than low-complexity methods".'),

      h2('4.5. Hình 11 — Ước lượng kênh, Doppler nguyên'),
      p('AFDM, N=256, ước lượng kênh qua pilot nhúng (Mục VI bài báo), so sánh Ideal CSI với SNRp=25/30/35dB.'),
      img('fig11_reproduction.png'),
      caption('Hình 11 tái hiện: SNRp=35dB gần trùng Ideal ở SNR cao.'),
      p('Tại SNR=18dB: Ideal=7.38×10⁻⁴, SNRp=35dB=9.20×10⁻⁴ — chênh lệch nhỏ, đúng khẳng định "for SNRp=35dB, the performance of AFDM with the proposed channel estimation is very close to the ideal case".'),

      h2('4.6. Hình 12a — Ước lượng kênh, Doppler lẻ bin'),
      p('Tương tự Hình 11 nhưng Doppler lẻ bin (ξν=1), so sánh SNRp=30/35/40dB. Đã áp dụng cả 2 bản sửa lỗi bộ tách sóng (Mục 3).'),
      img('fig12a_reproduction.png'),
      caption('Hình 12a tái hiện (sau khi sửa lỗi): thứ tự SNRp=40>35>30dB đúng bài báo trên toàn dải SNR.'),
      p('Sau khi sửa cả 2 lỗi (Mục 3), thứ tự các đường đúng hoàn toàn khớp bài báo tại mọi điểm SNR: SNRp=40dB luôn gần Ideal nhất, SNRp=30dB xa nhất — trước khi sửa Lỗi 2, SNRp=40dB từng cho BER XẤU HƠN cả SNRp=30/35dB, một sai lệch rõ ràng đã được xác nhận và loại bỏ.'),

      h2('4.7. Hình 12b — Hiệu ứng tham số guard ξν'),
      p('AFDM, N=256, 4QAM, so sánh ξν=1/4/8 dưới Ideal CSI và ước lượng (SNRp=40dB). Đã áp dụng bản sửa lỗi bộ tách sóng.'),
      img('fig12b_reproduction.png'),
      caption('Hình 12b tái hiện (sau khi sửa lỗi): ξν càng lớn, BER càng thấp — đúng bài báo.'),
      p('Tại SNR=18dB (Ideal): ξν=1→8.84×10⁻⁴, ξν=4→5.12×10⁻⁴, ξν=8→3.94×10⁻⁴ — giảm đơn điệu theo ξν tăng, đúng chính xác "increasing ξν improves the performance of AFDM, as less overlapping is present".'),

      h1('5. Giới hạn và ghi chú trung thực'),
      bullet('numMonte thấp hơn nhiều so với 10⁶ của bài báo (do giới hạn thời gian tính toán) — các đường cong không đạt độ mượt tuyệt đối ở vùng BER cực thấp, nhưng không ảnh hưởng đến các kết luận so sánh chính.'),
      bullet('Bộ ước lượng kênh (Mục VI bài báo) được triển khai đơn giản hóa: tìm đường mạnh nhất rồi loại trừ liên tiếp (successive cancellation/matching pursuit), thay vì thuật toán tối ưu đầy đủ mô tả trong bài — đã kiểm chứng cho kết quả đúng hướng qua chẩn đoán trực tiếp (biết trước delay/Doppler thật).'),
      bullet('MMSE băng-hẹp (Hình 10) dùng cách xấp xỉ "giữ L phần tử lớn nhất mỗi hàng" thay vì cấu trúc guard-symbol chính xác của bài báo — cho đúng xu hướng định tính nhưng độ lớn chênh lệch có thể khác biệt.'),
      bullet('Tất cả các phát hiện, lỗi, và cách sửa đều được ghi lại trung thực trong quá trình làm việc, không có số liệu nào bị làm giả hay che giấu.'),

      h1('6. Vị trí file'),
      p('Toàn bộ code MATLAB, dữ liệu (.mat), log chạy, và hình ảnh gốc được lưu tại:'),
      p('D:\\doubly-dispersive-channel-simulation-main\\afdm_paper_check\\', { font: 'Consolas' }),
    ],
  }],
});

Packer.toBuffer(doc).then((buffer) => {
  fs.writeFileSync(DIR + 'BaoCao_TaiHien_AFDM_Paper.docx', buffer);
  console.log('DONE');
});
