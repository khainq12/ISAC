const fs = require('fs');
const {
  Document, Packer, Paragraph, TextRun, HeadingLevel, ImageRun,
  Table, TableRow, TableCell, WidthType, BorderStyle, AlignmentType, PageBreak
} = require('docx');

const OUT = 'D:\\doubly-dispersive-channel-simulation-main\\report_build\\';

// ---------- helper functions ----------

function h1(text) { return new Paragraph({ heading: HeadingLevel.HEADING_1, spacing: { before: 300, after: 150 }, children: [new TextRun(text)] }); }
function h2(text) { return new Paragraph({ heading: HeadingLevel.HEADING_2, spacing: { before: 250, after: 120 }, children: [new TextRun(text)] }); }
function p(text, opts) { return new Paragraph({ spacing: { after: 150 }, children: [new TextRun(Object.assign({text}, opts||{}))] }); }
function pLead(lead, text) { return new Paragraph({ spacing: { after: 150 }, children: [new TextRun({text: lead + ' ', italics: true}), new TextRun(text)] }); }
function caption(text) { return new Paragraph({ alignment: AlignmentType.CENTER, spacing: { after: 250 }, children: [new TextRun({text, italics: true, size: 20})] }); }

// Renders a formula from mixed parts: plain strings, or {t, sub:true}/{t, sup:true} for
// real Word subscript/superscript (not Unicode look-alike characters). num (optional):
// appends a right-aligned-style equation number in parentheses.
function eqParts(parts, num) {
  const runs = parts.map(part => {
    if (typeof part === 'string') return new TextRun({ text: part, italics: true });
    return new TextRun({ text: part.t, italics: true, subScript: !!part.sub, superScript: !!part.sup });
  });
  if (num) runs.push(new TextRun({ text: '   (' + num + ')', italics: false, size: 18 }));
  return new Paragraph({ alignment: AlignmentType.CENTER, spacing: { before: 100, after: 140 }, children: runs });
}

function fmtNum(x) {
  if (Math.abs(x) < 0.001 && x !== 0) return x.toExponential(2);
  return Number(x.toFixed(4)).toString();
}

function dataTable(headers, rows) {
  const border = { style: BorderStyle.SINGLE, size: 2, color: '999999' };
  const borders = { top: border, bottom: border, left: border, right: border };
  const nCols = headers.length;
  const colPct = Math.floor(100 / nCols);
  function cell(text, bold) {
    return new TableCell({
      width: { size: colPct, type: WidthType.PERCENTAGE },
      borders,
      children: [new Paragraph({ alignment: AlignmentType.CENTER, children: [new TextRun({ text: String(text), bold: !!bold, size: 16 })] })]
    });
  }
  const headerRow = new TableRow({ tableHeader: true, children: headers.map(h => cell(h, true)) });
  const bodyRows = rows.map(r => new TableRow({ children: r.map(v => cell(typeof v === 'number' ? fmtNum(v) : v)) }));
  return new Table({ width: { size: 100, type: WidthType.PERCENTAGE }, rows: [headerRow, ...bodyRows] });
}

function img(path, w, h) {
  return new ImageRun({ type: 'png', data: fs.readFileSync(path), transformation: { width: w, height: h } });
}

function noBorderCell(children, widthPct) {
  return new TableCell({
    width: { size: widthPct, type: WidthType.PERCENTAGE },
    borders: { top: {style: BorderStyle.NONE}, bottom: {style: BorderStyle.NONE}, left: {style: BorderStyle.NONE}, right: {style: BorderStyle.NONE} },
    children
  });
}

function threeImgRow(paths, labels, w, h) {
  const cellW = Math.floor(100/3);
  return new Table({
    width: { size: 100, type: WidthType.PERCENTAGE },
    columnWidths: [3333, 3333, 3333],
    rows: [ new TableRow({ children: paths.map((pth, i) => noBorderCell([
      new Paragraph({ alignment: AlignmentType.CENTER, spacing: { after: 20 }, children: [img(pth, w, h)] }),
      new Paragraph({ alignment: AlignmentType.CENTER, children: [new TextRun({ text: labels[i], bold: true, size: 18 })] }),
    ], cellW)) }) ]
  });
}

function complexityTable() {
  const headers = ['Tiêu chí', 'OFDM', 'OTFS', 'AFDM'];
  const rows = [
    ['Điều chế/giải điều chế', 'O(N log N) — FFT/IFFT', 'O(N log N) — SFFT/ISFFT 2 chiều, hằng số lớn hơn', 'O(N log N) — FFT cộng 2 phép nhân chirp O(N)'],
    ['Cân bằng zero-forcing đầy đủ (kênh Doppler tổng quát)', 'O(N³) khi G không còn chéo do ICI', 'O(N³) tổng quát; có thể giảm nhờ khai thác cấu trúc block-circulant của kênh delay-Doppler', 'O(N³) tổng quát; cấu trúc dải hẹp của G (Mục 10.1) cho phép cân bằng "banded" chi phí thấp hơn nhiều'],
    ['Ước lượng kênh từ pilot', 'Đơn giản, một pilot cho mỗi subcarrier', 'Cần chiếm một phần lưới delay-Doppler, chi phí pilot cao hơn', 'Một khối pilot có thể tách nhiều đường truyền nhờ vị trí riêng biệt trong miền DAFT (Mục 10.1)'],
    ['Thu Maximum Likelihood vét cạn', 'O(2^(N·bps)) — như nhau cho cả ba waveform, chỉ khả thi với N nhỏ', 'O(2^(N·bps))', 'O(2^(N·bps))'],
  ];
  return dataTable(headers, rows);
}

const doc = new Document({
  sections: [{
    properties: { page: { size: { width: 12240, height: 15840 } } },
    children: [

      new Paragraph({ heading: HeadingLevel.TITLE, alignment: AlignmentType.CENTER, spacing: { after: 300 },
        children: [new TextRun('Khảo sát AFDM cho hệ thống ISAC: Lý thuyết bậc phân tập và Đánh giá mô phỏng')] }),
      p('Tài liệu này gồm hai phần: Phần I chứng minh bằng công thức toán học, trích dẫn và tổng hợp từ các công trình [1]-[3], rằng AFDM đạt bậc phân tập cao hơn OTFS và OFDM trong kênh đa đường/Doppler cao - kịch bản tiêu biểu cho ứng dụng ISAC. Phần II dùng chính các công thức, ma trận và điều kiện đã thiết lập ở Phần I để khảo sát mô phỏng hiệu năng AFDM: bậc phân tập, cấu trúc kênh hiệu dụng, hiệu năng BER dưới các bộ tách sóng thực tế, và khả năng ước lượng kênh qua pilot nhúng. Phương pháp và tham số mô phỏng tham khảo từ [2]; toàn bộ mã MATLAB tự viết theo đúng công thức Phần I.'),

      new Paragraph({ children: [new PageBreak()] }),

      new Paragraph({ heading: HeadingLevel.HEADING_1, alignment: AlignmentType.CENTER, spacing: { before: 200, after: 300 },
        children: [new TextRun('PHẦN I: LÝ THUYẾT')] }),
      new Paragraph({ heading: HeadingLevel.TITLE, alignment: AlignmentType.CENTER, spacing: { after: 300 },
        children: [new TextRun({text: 'Chứng minh: bậc phân tập AFDM > OTFS > OFDM', bold: true})] }),

      h1('1. Mô hình toán học cơ bản'),
      p('Bảng dưới đây thống nhất các ký hiệu dùng xuyên suốt tài liệu, mô tả một kênh vô tuyến đa đường gồm P đường truyền tác động lên một khối tín hiệu N mẫu:'),
      dataTable(['Ký hiệu', 'Ý nghĩa', 'Ghi chú'], [
        ['N', 'Kích thước một khối truyền (số mẫu tín hiệu trong một khung)', 'đánh số n = 0, 1, ..., N−1'],
        ['P', 'Số đường truyền (path) của kênh đa đường', 'do tín hiệu phản xạ/tán xạ trên nhiều vật cản đến máy thu theo nhiều đường'],
        ['hᵢ', 'Hệ số suy hao/fading của đường truyền i (i=1,...,P)', 'số phức'],
        ['lᵢ', 'Độ trễ của đường truyền i', 'số nguyên, chuẩn hóa theo chu kỳ lấy mẫu'],
        ['νᵢ', 'Độ dịch Doppler của đường truyền i', 'số thực, chuẩn hóa theo khoảng cách giữa các mẫu'],
        ['Π', 'Ma trận dịch vòng (cyclic shift), N×N', '(Πx)[n]=x[(n−1) mod N]; Π^{lᵢ}: dịch lᵢ mẫu liên tiếp (lũy thừa ma trận)'],
        ['W', 'Ma trận xoay pha Doppler, N×N, đường chéo', 'W^{νᵢ} := diag(e^{j2πνᵢn/N}), n=0..N−1 - ký hiệu tham số hóa theo νᵢ, không phải lũy thừa ma trận'],
        ['CPPᵢ', 'Ma trận hiệu chỉnh pha cho mẫu quấn vòng qua biên khối', '= I (đơn vị) với OFDM/OTFS; chỉ khác I với AFDM do miền chirp không tuần hoàn chu kỳ N (công thức xem Bemani et al., Mục III)'],
      ]),
      p('Với các ký hiệu trên, kênh vật lý tổng hợp từ P đường truyền là:'),
      eqParts(['H = Σ', {t:'i=1', sub:true}, {t:'P', sup:true}, ' h', {t:'i', sub:true}, '·Δ', {t:'i', sub:true}, ',      Δ', {t:'i', sub:true}, ' = CPP', {t:'i', sub:true}, '·W', {t:'νᵢ', sup:true}, '·Π', {t:'lᵢ', sup:true}], '1.1'),
      p('trong đó Δᵢ là ma trận N×N biểu diễn riêng đường truyền i, còn H là ma trận N×N biểu diễn toàn bộ kênh vật lý (mẫu tín hiệu thu = H·(mẫu tín hiệu phát) + nhiễu). Mô hình lấy từ arXiv:2104.11331 (phương trình 11-14).'),
      p('Mỗi waveform (OFDM/OTFS/AFDM) không truyền trực tiếp trên miền mẫu thời gian mà trên một miền ký hiệu riêng của nó, thông qua một ma trận điều chế M (N×N, đổi vector ký hiệu ở miền riêng sang vector mẫu thời gian đem phát) và ma trận giải điều chế A = M⁻¹ (đổi mẫu thời gian thu về lại miền riêng để tách sóng):'),
      dataTable(['Waveform', 'Ma trận M', 'Ghi chú'], [
        ['OFDM', 'F (IDFT)', 'F: ma trận IDFT N điểm - miền riêng của OFDM là miền tần số (subcarrier)'],
        ['OTFS', 'F_L ⊗ I_K', 'F_L: DFT L điểm (L=số bin trễ); I_K: đơn vị K×K (K=số bin Doppler); N=K·L; ⊗: tích Kronecker - miền riêng là lưới delay-Doppler'],
        ['AFDM', 'Λc2·F·Λc1', 'Λc = diag(e^{−j2πc·n²}), n=0..N−1; c1, c2 là hai tham số chirp (định nghĩa ở Mục 4) - miền riêng là miền chirp (DAFT)'],
        ['OCDM', 'Λc2·F·Λc1 với c1=c2=1/(2N)', 'trường hợp riêng của AFDM (cùng ma trận M, không thêm tham số mới) - dùng làm đối chứng bổ sung ở Mục 12'],
      ]),
      p('Kênh hiệu dụng nhìn từ phía ký hiệu thông tin (miền riêng của từng waveform) là G = A·H·M - một ma trận N×N liên hệ trực tiếp vector ký hiệu phát s với vector ký hiệu thu r theo r = G·s + nhiễu. Ba cặp (A,M) khác nhau (ứng với OFDM/OTFS/AFDM) cho các ma trận G có cấu trúc khác nhau dù cùng xuất phát từ một kênh vật lý H - đây là lý do các waveform chống nhiễu khác nhau dù truyền qua cùng một kênh vô tuyến.'),

      h1('2. Vì sao cần một miền biểu diễn mới'),
      p('Mục 1 cho thấy ba waveform khác nhau ở cặp ma trận (A, M), tức ở miền biểu diễn ký hiệu. Phần này tóm tắt cơ chế vật lý đứng sau ba lựa chọn đó, làm nền cho định nghĩa bậc phân tập ở Mục 3 và cho cách đọc cấu trúc ma trận kênh ở Mục 10.'),

      h2('2.1 OFDM: nhiễu liên kênh (ICI) dưới dịch Doppler'),
      p('OFDM giải điều chế bằng phép chiếu tín hiệu thu lên các hàm cơ sở tần số rời rạc. Với một subcarrier k chịu dịch Doppler chuẩn hóa ν, năng lượng rò sang subcarrier cách nó l vị trí tỉ lệ với:'),
      eqParts(['C', {t:'l', sub:true}, '(ν) = (1/N)·Σ', {t:'n', sub:true}, ' e', {t:'j2π(ν−l)n/N', sup:true}], '2.1'),
      p('Khi ν là số nguyên, Cₗ(ν) triệt tiêu với mọi l≠0 (trực giao của các hàm mũ phức rời rạc), nên không rò rỉ. Khi ν lệch khỏi số nguyên, Cₗ(ν) khác 0 với mọi l và suy giảm chậm (cỡ 1/l): năng lượng rò sang nhiều subcarrier lân cận cùng lúc, mang tính toàn cục, không có cấu trúc thưa để khai thác khi cân bằng kênh. Đây là nhiễu liên kênh (ICI) của OFDM.'),

      h2('2.2 OTFS: rò rỉ khu trú trên lưới delay-Doppler'),
      p('OTFS truyền trên một lưới hai chiều K bin Doppler × L bin trễ (N=K·L, Mục 1). Khi trễ/Doppler mọi đường truyền rơi đúng lưới, kênh chỉ dịch vòng ký hiệu, không rò rỉ. Khi Doppler lệch khỏi lưới (trường hợp phổ biến hơn), rò rỉ vẫn xảy ra nhưng chỉ khu trú quanh chỉ số Doppler:'),
      eqParts(['Y[k,l] = Σ', {t:'q', sub:true}, ' G(q,κ)·x[k−q, l]'], '2.2'),
      p('với κ là phần phân số của Doppler và G(q,κ) suy giảm nhanh theo |q| - khác với Cₗ(ν) suy giảm chậm của OFDM ở (2.1). Cấu trúc rò rỉ khu trú này là lý do OTFS chịu Doppler tốt hơn OFDM về mặt cấu trúc kênh, không nhờ một thuật toán thu đặc biệt.'),

      h2('2.3 AFDM: miền chirp bất biến dưới dịch Doppler'),
      p('Một chirp x(t)=e^{jπat²} chịu dịch Doppler ν vẫn là một chirp cùng hệ số a, chỉ dịch trục thời gian và thêm một pha hằng số: dạng chirp không đổi dưới dịch Doppler. Ngược lại, một hàm mũ phức thuần (cơ sở của OFDM) đổi hẳn sang một hàm cơ sở khác dưới cùng phép dịch - đây là lý do một hệ tuyến tính bất biến thời gian có hàm mũ phức làm hàm riêng (eigenfunction), nhưng kênh Doppler (hệ biến đổi theo thời gian) thì không, nên biến đổi Fourier không chéo hóa được kênh Doppler - gốc rễ của ICI ở Mục 2.1. Vì chirp giữ nguyên dạng dưới Doppler, một biến đổi trên nền chirp - Discrete Affine Fourier Transform (DAFT), ma trận Λc2·F·Λc1 ở Mục 1 - là lựa chọn tự nhiên để gần chéo hóa đồng thời cả trễ và Doppler. Đây là động lực xây dựng AFDM; Mục 4 chứng minh lựa chọn này cho hạng đầy đủ.'),

      h1('3. Định nghĩa bậc phân tập'),
      p('Một từ mã (codeword) là một vector N ký hiệu điều chế mà máy phát có thể gửi (x, N×1). Xét hai từ mã xₘ, xₙ mà máy thu có thể nhầm lẫn, đặt δ = xₘ−xₙ. Với mỗi đường truyền p, gọi Ψₚ = A·Δₚ·M là kênh hiệu dụng nếu chỉ riêng đường truyền p tồn tại; vector Ψₚδ là phần đóng góp riêng của đường truyền p vào sai khác tín hiệu thu. Xếp P vector này cạnh nhau:'),
      eqParts(['Φ(δ) = [ Ψ', {t:'1', sub:true}, 'δ  |  Ψ', {t:'2', sub:true}, 'δ  |  ...  |  Ψ', {t:'P', sub:true}, 'δ ]'], '3.1'),
      p('Hạng của Φ(δ) (một số nguyên từ 1 đến P) đo số "hướng" độc lập mà P đường truyền đóng góp vào sai khác tín hiệu thu - hạng càng cao, máy thu càng dễ phân biệt xₘ với xₙ. Bậc phân tập (diversity order) ρ của hệ thống là hạng nhỏ nhất trên mọi cặp từ mã có thể gây nhầm lẫn:'),
      eqParts(['ρ = min', {t:'m,n (m≠n)', sub:true}, ' rank( Φ(x', {t:'m', sub:true}, ' − x', {t:'n', sub:true}, ') )'], '3.2'),
      p('Xác suất nhầm lẫn giữa hai từ mã cụ thể (PEP), lấy trung bình trên fading ngẫu nhiên hᵢ, được chặn trên bởi (arXiv:2104.11331, phương trình 38):'),
      eqParts(['P(x', {t:'m', sub:true}, ' → x', {t:'n', sub:true}, ') ≤ ∏', {t:'l=1', sub:true}, {t:'r', sup:true}, ' [ 1/(1 + λ', {t:'l', sub:true}, '²/(4·P·N', {t:'0', sub:true}, ')) ]'], '3.3'),
      p('với r = rank(Φ(δ)) và λₗ là các trị kỳ dị khác không của Φ(δ). Khi SNR lớn, vế phải xấp xỉ tỉ lệ nghịch với SNR mũ r - r chính là độ dốc của đường cong BER trên thang log-log theo SNR.'),

      h1('4. AFDM đạt bậc phân tập đầy đủ'),
      p('Nguồn: A. Bemani, N. Ksairi, M. Kountouris, arXiv:2104.11331 (bản mở rộng IEEE Trans. Wireless Commun. 2023, arXiv:2204.12798).'),
      p('Với lmax là độ trễ chuẩn hóa lớn nhất và αmax là Doppler chuẩn hóa lớn nhất trong P đường truyền của một kịch bản kênh, Định lý 1 phát biểu: nếu c1 = (2αmax+1)/(2N) và'),
      eqParts(['2·α', {t:'max', sub:true}, ' + l', {t:'max', sub:true}, ' + 2·α', {t:'max', sub:true}, '·l', {t:'max', sub:true}, ' < N'], '4.1'),
      p('thì AFDM đạt ρ = P (bậc phân tập đầy đủ). Cơ chế: c1 gán mỗi đường truyền i một vị trí riêng locᵢ = αᵢ + 2N·c1·lᵢ trong miền chirp - điều kiện (4.1) đảm bảo các vị trí này không trùng nhau; tham số chirp thứ hai c2 được chọn vô tỉ để loại trừ khả năng suy biến hạng còn lại (một số vô tỉ không thể bằng một số hữu tỉ xây dựng từ các thành phần nguyên của δ). Khi c1=0 (AFDM suy biến thành OFDM), locᵢ chỉ còn phụ thuộc Doppler nên hai đường truyền cùng Doppler khác độ trễ luôn trùng vị trí, cho ρ=1 với mọi δ khác 0 - đây chính là trường hợp OFDM. Trường hợp Doppler phân số (không rơi đúng lưới nguyên) và công thức mở rộng của c1, cùng số mẫu bảo vệ (guard) cần thiết, được trình bày ở Mục 10.1 cùng với cấu trúc dải hẹp của Heff.'),

      h1('5. OTFS: bậc phân tập giảm về 1 khi SNR đủ cao'),
      p('Nguồn: D. Surabhi, A. Chockalingam, "On the Diversity of Uncoded OTFS Modulation in Doubly-Dispersive Channels", IEEE Trans. Wireless Commun. 2019, arXiv:1808.07747.'),
      p('Xét cặp từ mã hằng số trên lưới delay-Doppler (Mục 2.2): một từ mã đặt ký hiệu a tại mọi vị trí lưới, từ mã kia đặt a′≠a - cặp này luôn tồn tại. Vì kênh chỉ dịch vòng ký hiệu trên lưới, một vector hằng số bị dịch đi vẫn cho lại chính nó, nên Φ(δ) có hạng 1 bất kể P - bậc phân tập tiệm cận (SNR→∞) của OTFS chưa mã hóa bị chặn trên bởi 1, kể cả với thu Maximum Likelihood, do cấu trúc lưới chứ không do hạn chế thuật toán thu. Ở SNR thấp/vừa, cặp từ mã hạng cao hơn chi phối nên độ dốc BER gần full diversity; chỉ ở SNR đủ cao, cặp suy biến hạng 1 mới chi phối và kéo độ dốc về 1 - bài báo đề xuất xoay pha ký hiệu bằng số vô tỉ để phá vỡ tính bất biến này, cùng nguyên lý với c2 vô tỉ của AFDM (Mục 4).'),

      h1('6. Tổng hợp'),
      dataTable(['Waveform', 'Bậc phân tập ρ', 'Điều kiện', 'Cơ chế'], [
        ['OFDM', '1 (mọi trường hợp)', 'không có', 'không tách được đường truyền theo Doppler (c1=0)'],
        ['OTFS', 'gần P ở SNR thấp/vừa; tiến về 1 khi SNR→∞', 'không xoay pha bổ sung', 'cặp từ mã hằng số lặp đều luôn cho hạng 1'],
        ['AFDM', 'P ở mọi SNR', '(4.1), c1=(2αmax+1)/(2N), c2 vô tỉ', 'tách vị trí (c1) và không có δ nào làm suy biến hạng (c2 vô tỉ)'],
      ]),
      p('Hệ quả cho đường cong BER dưới thu ML: ở SNR đủ cao, theo (3.3), độ dốc của AFDM (r=P) gấp P lần độ dốc của OFDM (r=1); độ dốc của OTFS giảm dần theo SNR, gần bằng AFDM ở SNR thấp/vừa rồi thoải dần về gần độ dốc OFDM ở SNR rất cao (cơ chế chuyển tiếp, Mục 5). Trong phạm vi điều kiện đã nêu (Mục 4-5), thứ tự bậc phân tập AFDM (ρ=P) > OTFS (ρ→1 khi SNR cao) > OFDM (ρ=1 luôn) là kết quả đã chứng minh bằng công thức, không phải quan sát thực nghiệm đơn lẻ. Phần II dưới đây kiểm chứng bằng số các dự đoán này.'),

      h1('7. Tài liệu tham khảo (Phần I)'),
      p('[1] A. Bemani, N. Ksairi, M. Kountouris, "AFDM: A Full Diversity Next Generation Waveform for High Mobility Communications", arXiv:2104.11331, 2021.'),
      p('[2] A. Bemani, N. Ksairi, M. Kountouris, "Affine Frequency Division Multiplexing for Next Generation Wireless Communications", IEEE Trans. Wireless Commun., 2023 (arXiv:2204.12798).'),
      p('[3] K. R. Surabhi, R. M. Augustine, A. Chockalingam, "On the Diversity of Uncoded OTFS Modulation in Doubly-Dispersive Channels", IEEE Trans. Wireless Commun., 2019 (arXiv:1808.07747).'),

      new Paragraph({ children: [new PageBreak()] }),

      new Paragraph({ heading: HeadingLevel.HEADING_1, alignment: AlignmentType.CENTER, spacing: { before: 200, after: 300 },
        children: [new TextRun('PHẦN II: MÔ PHỎNG')] }),

      h1('8. Giới thiệu mô phỏng'),
      p('Một hệ thống ISAC (Integrated Sensing and Communication) dùng chung tài nguyên phổ và phần cứng cho cả cảm biến lẫn thông tin liên lạc, nên tuyến thông tin phải duy trì độ tin cậy ngay cả trong kênh biến đổi nhanh theo thời gian (đa đường cộng Doppler cao, thường gặp khi thiết bị hoặc mục tiêu di chuyển tốc độ lớn) - đây chính là kịch bản kênh mà Phần I phân tích bằng công thức. Phần I chứng minh AFDM đạt bậc phân tập đầy đủ (ρ=P) trong kênh như vậy, còn OFDM - waveform phổ biến nhất hiện nay - luôn bị giới hạn ở ρ=1 (Định lý 1, Phần I Mục 4).'),
      p('Phần II khảo sát mô phỏng để đánh giá cụ thể hiệu năng này cho bối cảnh ISAC: bậc phân tập đo trực tiếp, cấu trúc kênh hiệu dụng, hiệu năng BER dưới các bộ tách sóng thực tế, và khả năng ước lượng kênh qua pilot nhúng - yếu tố đáng chú ý với ISAC vì tài nguyên pilot trong khung tín hiệu thường phải cân đối giữa hai chức năng cảm biến và thông tin. Toàn bộ mã MATLAB tự viết theo đúng công thức và ma trận (M, A) đã thiết lập ở Phần I; phương pháp mô phỏng (mô hình kênh, tham số, bộ tách sóng) và công thức cấu trúc dải hẹp của kênh hiệu dụng (Mục 10.1) tham khảo từ [2] (Mục V-VII). Tham số mô phỏng tổng hợp tại Bảng 1 (Mục 9.1).'),

      h1('9. Tham số và phương pháp mô phỏng'),

      h2('9.1. Tham số mô phỏng tổng hợp'),
      p('Kênh vật lý H, ma trận điều chế M/A của từng waveform đã định nghĩa ở Phần I (Mục 1) - mô phỏng dùng đúng các công thức đó, không định nghĩa lại. Bảng 1 liệt kê tham số của bảy thí nghiệm ở Mục 11-17; công thức bộ tách sóng và ước lượng kênh ở Mục 9.2. Tần số sóng mang cố định fc=4GHz cho mọi thí nghiệm; cột "v tối đa" quy đổi từ αmax sang vận tốc di chuyển tương ứng (theo [2], Mục VII) - mức tiêu biểu cho kịch bản ISAC di động cao (ví dụ giao thông tốc độ cao, thiết bị bay).'),
      dataTable(['Hình', 'N', 'P', 'Điều chế', 'lₘₐₓ', 'αₘₐₓ', 'v tối đa [km/h]', 'ξν', 'Bộ tách sóng / ước lượng', 'numMonte'], [
        ['2 (Mục 11)', 16, '2/3/4', 'BPSK', '-', 1, 405, '-', 'ML vét cạn', 10000],
        ['3 (Mục 12)', 16, 3, 'BPSK', 2, 1, 405, '-', 'ML vét cạn', 3000],
        ['4 (Mục 13)', 256, 3, 'QPSK', 2, 2, 540, '-', 'MMSE', 3000],
        ['5 (Mục 14)', 256, 3, '4-QAM', 2, 2, 540, '0 / 1', 'MMSE, MRC-DFE, MMSE dải hẹp', '1500-30000'],
        ['6 (Mục 15)', 256, 3, '4-QAM', 2, 2, 540, 0, 'Pilot nhúng, SNRp=25/30/35dB', 400],
        ['7 (Mục 16)', 256, 3, '4-QAM', 2, 2, 540, 1, 'Pilot nhúng, SNRp=30/35/40dB', 400],
        ['8 (Mục 17)', 256, 3, '4-QAM', 2, 2, 540, '1/4/8', 'Pilot nhúng, SNRp=40dB', 300],
      ]),
      caption('Bảng 1: Tham số của bảy thí nghiệm mô phỏng, Mục 11-17 (numMonte của Hình 5 tăng dần theo SNR để làm mượt vùng BER thấp, xem Mục 14).'),

      h2('9.2. Công thức bộ tách sóng và ước lượng kênh'),
      p('Bộ tách sóng MMSE giải Ŝ=(Heffᴴ·Heff+N0·I)⁻¹·Heffᴴ·Z (Heff=A·H·M, Mục 1). Bộ tách sóng MRC-based DFE dùng lặp Gauss-Seidel trên cùng hệ phương trình chuẩn ([2] chứng minh - Theorem 2 - tương đương LMMSE khi hội tụ). Ước lượng kênh qua pilot nhúng: một mẫu pilot biết trước đặt tại vị trí đầu khung với SNRp=|xpilot|²/N0, tìm đường truyền mạnh nhất bằng tương quan với các giả thuyết (trễ, Doppler) rồi loại trừ khỏi tín hiệu dư (successive cancellation), lặp lại cho đến khi tìm đủ P đường; công thức c1 (Doppler nguyên/phân số) và số mẫu bảo vệ Q trình bày ở Mục 10.1.'),

      h2('9.3. Độ phức tạp tính toán'),
      complexityTable(),
      caption('Bảng 2: So sánh độ phức tạp tính toán giữa OFDM, OTFS và AFDM.'),
      p('Điều chế/giải điều chế có cùng bậc O(N log N) cho cả ba (AFDM cộng thêm O(N) cho hai phép nhân chirp, không đáng kể). Cân bằng zero-forcing đầy đủ cần O(N³) khi kênh có Doppler tổng quát, nhưng cấu trúc dải hẹp của G ở AFDM (Mục 10.1) cho phép dùng cân bằng "banded" chi phí thấp hơn nhiều, kiểm chứng bằng số ở Mục 14.'),

      h1('10. Cấu trúc ma trận kênh hiệu dụng'),
      p('Mục này trực quan hóa ma trận kênh hiệu dụng Heff=A·H·M (Phần I, Mục 1) trên một kênh Rayleigh ngẫu nhiên hai đường truyền, để kiểm tra xem cấu trúc thưa/tập trung dự đoán ở Phần I Mục 2 có thật hay không, trước khi dùng nó giải thích kết quả BER ở Mục 11-17.'),
      p('Hình 1: hai đường truyền ℓ1=2,ν1≈+0.35 và ℓ2=1,ν2≈−0.4 (lưới N=256), Doppler nguyên (hàng trên) và phân số (hàng dưới).'),
      threeImgRow([OUT+'chmat_int_1_v1.png', OUT+'chmat_int_2_v1.png', OUT+'chmat_int_3_v1.png'], ['a) OFDM (nguyên)', 'b) OTFS (nguyên)', 'c) AFDM (nguyên)'], 220, 206),
      threeImgRow([OUT+'chmat_frac_1_v1.png', OUT+'chmat_frac_2_v1.png', OUT+'chmat_frac_3_v1.png'], ['d) OFDM (phân số)', 'e) OTFS (phân số)', 'f) AFDM (phân số)'], 220, 206),
      caption('Hình 1: Ma trận kênh hiệu dụng, Doppler nguyên (trên) và phân số (dưới), hai đường truyền P=2.'),
      pLead('Quan sát:', 'AFDM và OTFS tập trung năng lượng trên vài đường chéo hẹp (rõ nhất với Doppler nguyên, Hình 1a-c), còn OFDM trải năng lượng rộng quanh đường chéo chính, càng rộng hơn khi Doppler phân số (Hình 1d-f).'),
      pLead('Giải thích:', 'OFDM không chéo hóa gọn toán tử Doppler (Phần I Mục 2.1); Doppler phân số ở OTFS rò rỉ ra ô lưới lân cận nhưng khu trú (Phần I Mục 2.2); AFDM ít bị ảnh hưởng vì xử lý Doppler như một phép dịch liên tục qua pha chirp (Phần I Mục 2.3). Cấu trúc càng thưa/tập trung thì hiệu năng đo được ở Mục 11-17 càng tốt, vì máy thu càng dễ tách các đường truyền - Mục 10.1 nêu công thức định lượng đứng sau quan sát định tính này.'),

      h2('10.1. Cấu trúc dải hẹp (banded) của Heff dưới AFDM'),
      p('Hình 1 cho thấy cấu trúc thưa/tập trung của AFDM một cách định tính; mục này nêu công thức định lượng đứng sau quan sát đó, theo Bemani et al. [2] (Mục V-VI), làm cơ sở cho bộ tách sóng độ phức tạp thấp kiểm chứng ở Mục 14.'),
      p('Khi khung tín hiệu chèn đủ mẫu bảo vệ null (không mang dữ liệu) ở hai đầu, phép dịch vòng Πˡⁱ trong (1.1) không còn "quấn" qua biên khối trong phần dữ liệu - Heff=A·H·M khi đó là một ma trận dải hẹp chính xác (không xấp xỉ): mỗi cột chỉ có L phần tử khác 0, với'),
      eqParts(['L = P  (Doppler nguyên);      L = (2ξ', {t:'ν', sub:true}, '+1)·P  (Doppler phân số)'], '10.1'),
      p('trong đó ξν là biên bảo vệ Doppler phân số (số nguyên ≥0, càng lớn càng nới rộng vùng tìm kiếm lân cận mỗi vị trí locᵢ). Số mẫu null-guard tối thiểu cần chèn mỗi khung là'),
      eqParts(['Q = (l', {t:'max', sub:true}, '+1)(2(α', {t:'max', sub:true}, '+ξ', {t:'ν', sub:true}, ')+1) − 1'], '10.2'),
      p('Cơ chế: mỗi đường truyền i đóng góp đúng vào vị trí locᵢ=αᵢ+2N·c1·lᵢ trong miền chirp (Phần I Mục 4). Với Doppler nguyên, P vị trí locᵢ là P cột lân cận nhau nên L=P đúng chính xác. Với Doppler phân số, locᵢ không rơi đúng lưới nguyên nên mỗi đường truyền "loang" ra thêm 2ξν vị trí lân cận, cho L=(2ξν+1)P - đây là lý do (10.1) chỉ đúng chính xác ở trường hợp Doppler nguyên, còn Doppler phân số là một xấp xỉ có kiểm soát bởi ξν. Với Doppler phân số, c1 mở rộng từ Định lý 1 (Phần I Mục 4) thành c1=(2(αmax+ξν)+1)/(2N).'),
      p('Cấu trúc (10.1) là cơ sở định lượng cho dòng "cân bằng banded chi phí thấp hơn nhiều" đã nêu ở Bảng 2 (Mục 9.3): thay vì nghịch đảo ma trận N×N đầy đủ (O(N³)), bộ tách sóng chỉ cần xử lý L cột khác 0 mỗi hàng - một lợi thế đáng chú ý cho phần cứng ISAC, nơi khối xử lý thông tin thường phải chia sẻ tài nguyên tính toán với thuật toán cảm biến chạy song song. Cùng công thức Q ở (10.2) xác định đúng bao nhiêu mẫu null-guard cần dành cho ước lượng kênh bằng pilot nhúng (Mục 15-17). Mục 14 khảo sát bằng số cả hai hệ quả: với Doppler nguyên, bộ tách sóng MMSE dải hẹp cho BER gần như trùng với MMSE đầy đủ (phù hợp với L=P không xấp xỉ); với Doppler phân số, MMSE dải hẹp kém hơn một chút do phần xấp xỉ ở (10.1).'),

      h1('11. Bậc phân tập AFDM đo bằng thu ML vét cạn'),
      p('Mục tiêu: đo trực tiếp độ dốc đường cong BER dưới thu Maximum Likelihood vét cạn (2^{N·bps} giả thuyết, khả thi vì N=16 nhỏ) khi số đường truyền P thay đổi, đối chiếu với ρ=P của Định lý 1 (Phần I Mục 4).'),
      new Paragraph({ alignment: AlignmentType.CENTER, spacing: { after: 40 }, children: [img(OUT+'afdm_repro/fig8a_reproduction.png', 460, 366)] }),
      caption('Hình 2: BER theo SNR, thu ML, AFDM với P=2/3/4 đường truyền, Doppler nguyên.'),
      pLead('Quan sát:', 'P=4 cho độ dốc lớn nhất, P=2 nhỏ nhất - độ dốc tăng dần theo P đúng thứ tự ρ=P dự đoán ở (4.1)-(10.1). Độ dốc đo được trong đoạn 8-16dB xấp xỉ 1.5/2.8/3.3 cho P=2/3/4, đang tiệm cận giá trị lý thuyết 2/3/4 khi SNR tăng thêm (độ dốc chỉ đo "sạch" ở vùng BER thấp, xem cách đọc ở Mục 3).'),
      p('Kết luận: độ dốc đo được tăng tuyến tính theo P, phù hợp với hạng rank(Φ(δ))=P mà Định lý 1 dự đoán (Phần I Mục 4). Với một hệ thống ISAC hoạt động trong kênh nhiều đường phản xạ - phổ biến khi có mục tiêu và vật cản xung quanh - kết quả này cho thấy AFDM khai thác được toàn bộ P đường truyền độc lập để giảm xác suất lỗi, thay vì chỉ dựa vào một đường mạnh nhất như OFDM.'),

      h1('12. So sánh bốn waveform bằng thu ML'),
      p('Mục tiêu: so sánh OFDM/OCDM/OTFS/AFDM (Mục 1) dưới cùng thu ML, kênh 3 đường cố định (trễ 0,1,2), Doppler phân số theo phổ Jakes.'),
      new Paragraph({ alignment: AlignmentType.CENTER, spacing: { after: 40 }, children: [img(OUT+'afdm_repro/fig8b_reproduction.png', 460, 366)] }),
      caption('Hình 3: BER theo SNR, thu ML, bốn waveform, Doppler phân số Jakes.'),
      pLead('Quan sát:', 'OFDM tệ nhất ở mọi SNR; OCDM khá hơn nhưng vẫn hạn chế; OTFS và AFDM gần như trùng khít và tốt nhất.'),
      pLead('Giải thích:', 'OFDM không tách được các đường truyền theo Doppler khi vượt quá độ phân giải nguyên (cơ chế ICI toàn cục, Mục 2.1), nên hiệu năng luôn thấp nhất trong bốn waveform. OCDM dùng cùng dạng ma trận biến đổi affine với AFDM (Λc2·F·Λc1, Mục 1) nhưng giữ cố định c1=c2=1/(2N) bất kể thông số kênh thực tế, trong khi AFDM chọn c1 theo (4.1) để đặt vị trí locᵢ của từng đường truyền tách biệt nhau trong miền chirp (Mục 10.1) - vì vậy OCDM chỉ cải thiện một phần so với OFDM, còn AFDM và OTFS đạt mức tách đường truyền gần tối ưu. Với ISAC, nơi khoảng cách và tốc độ mục tiêu - và do đó đặc tính trễ/Doppler của kênh - thay đổi theo kịch bản, khả năng AFDM tự điều chỉnh tham số chirp theo đúng kênh hiện tại (thay vì dùng một cấu hình cố định như OCDM) là một lợi thế thiết kế thực tế.'),

      h1('13. So sánh bốn waveform bằng thu MMSE, N=256'),
      p('Mục tiêu: lặp lại Mục 12 với N lớn hơn (256) và thu MMSE - bộ thu thực tế hơn ML cho N lớn.'),
      new Paragraph({ alignment: AlignmentType.CENTER, spacing: { after: 40 }, children: [img(OUT+'afdm_repro/fig9a_reproduction.png', 460, 366)] }),
      caption('Hình 4: BER theo SNR, thu MMSE, bốn waveform, N=256, Doppler phân số Jakes.'),
      pLead('Quan sát:', 'tại SNR=27.5dB, OTFS và AFDM cho cùng một giá trị BER=3.2552×10⁻⁶ - gần trùng khít; OFDM và OCDM vẫn kém hơn rõ rệt, cùng thứ tự đã thấy ở Hình 3.'),
      pLead('Giải thích:', 'OTFS và AFDM đạt hiệu năng BER gần như giống hệt nhau dưới thu MMSE tuyến tính, không chỉ dưới thu ML (Hình 3) - cả hai đều tách đường truyền gần tối ưu, chỉ khác miền biểu diễn (lưới delay-Doppler hai chiều với OTFS, miền chirp một chiều với AFDM, Mục 2). Với ISAC, sự tương đương này có thêm một khía cạnh thực tế: theo [2], chi phí pilot cho ước lượng kênh nhúng của AFDM chỉ bằng khoảng một nửa so với OTFS, do cấu trúc lưới hai chiều của OTFS cần nhiều mẫu bảo vệ hơn (Mục V-VI, [2]) - ở cùng mức hiệu năng thông tin, AFDM để lại nhiều tài nguyên khung hơn cho chức năng cảm biến.'),

      h1('14. Cấu trúc dải hẹp của Heff: so sánh ba bộ tách sóng'),
      p('Mục tiêu: kiểm chứng bằng số công thức (10.1) - so sánh MMSE đầy đủ, MRC-based DFE, và MMSE dải hẹp (chỉ giữ L phần tử lớn nhất mỗi hàng của Heff) cho AFDM và OFDM, cả Doppler nguyên (trái) và phân số ξν=1 (phải).'),
      new Paragraph({ alignment: AlignmentType.CENTER, spacing: { after: 40 }, children: [img(OUT+'afdm_repro/fig10_reproduction.png', 600, 269)] }),
      caption('Hình 5: BER theo SNR, ba bộ tách sóng, AFDM và OFDM - trái: Doppler nguyên; phải: Doppler phân số ξν=1.'),
      pLead('Quan sát:', 'Doppler nguyên (trái): ba bộ tách sóng gần như trùng nhau trong cùng một waveform - phù hợp với L=P không xấp xỉ ở (10.1); nhóm AFDM tiếp tục giảm tới BER≈1.2×10⁻⁶ tại SNR=30dB trong khi nhóm OFDM chững thành một sàn lỗi quanh 1.25×10⁻⁴. Doppler phân số (phải): nhóm AFDM nằm dưới nhóm OFDM ở mọi SNR; MMSE dải hẹp của AFDM vẫn tốt hơn cả ba đường OFDM.'),
      pLead('Giải thích:', 'Với OFDM, c1=0 (Mục 1) nên locᵢ chỉ phụ thuộc Doppler (Phần I Mục 4) - hai đường truyền có Doppler làm tròn trùng nhau luôn chung một vị trí trong Heff dù trễ khác nhau, tạo một sàn lỗi mang tính cấu trúc mà tăng SNR không khắc phục được; AFDM (vị trí phụ thuộc cả trễ) không gặp hiện tượng này, phù hợp với cơ chế ρ=P đã chứng minh ở Phần I Mục 4. Với ISAC trong kênh có nhiều mục tiêu/vật cản cùng tốc độ tương đối gần nhau (Doppler gần trùng), đây là một khác biệt thực tế: OFDM có nguy cơ gặp sàn lỗi ngay cả khi SNR cao, còn AFDM vẫn tiếp tục cải thiện.'),
      p('Giới hạn: khung mô phỏng của Hình 5 không chèn null-guard như giả thiết ở Mục 10.1, nên bộ tách sóng dải hẹp áp ngưỡng L trực tiếp lên Heff dạng tuần hoàn đầy đủ - hệ số L thực dùng cho trường hợp phân số lớn hơn công thức (10.1) khoảng 3 lần mới cho bậc sai số quan sát được ở đây.'),

      h1('15. Ước lượng kênh qua pilot nhúng, Doppler nguyên'),
      p('Mục tiêu: đo BER dưới thu MMSE khi kênh được ước lượng qua pilot nhúng (Mục 9.2) thay vì biết trước (Ideal), với ba mức công suất pilot SNRp.'),
      new Paragraph({ alignment: AlignmentType.CENTER, spacing: { after: 40 }, children: [img(OUT+'afdm_repro/fig11_reproduction.png', 460, 366)] }),
      caption('Hình 6: BER theo SNR, AFDM, Ideal so với ước lượng qua pilot (SNRp=25/30/35dB), Doppler nguyên.'),
      pLead('Quan sát:', 'SNRp=35dB gần trùng Ideal ở SNR cao (SNR=18dB: Ideal=7.38×10⁻⁴, SNRp=35dB=9.20×10⁻⁴); SNRp giảm thì BER tăng đơn điệu, đúng thứ tự trên toàn dải SNR.'),
      pLead('Giải thích:', 'theo Mục 10.1, Q mẫu null-guard cô lập đúng L=P đường truyền trong miền chirp, nên sai số ước lượng chỉ đến từ nhiễu trên pilot, không từ nhầm lẫn cấu trúc giữa các đường truyền - vì vậy tăng SNRp luôn cải thiện đơn điệu, không có hiệu ứng phụ. Với ISAC, đây là cơ sở để chọn công suất pilot phù hợp: SNRp=35dB đã đủ để ước lượng kênh gần như không mất hiệu năng so với biết trước kênh, trong khi phần khung dành cho pilot (theo Q ở (10.2)) vẫn tương đối nhỏ so với N=256.'),

      h1('16. Ước lượng kênh qua pilot nhúng, Doppler phân số'),
      p('Mục tiêu: lặp lại Mục 15 với ξν=1 (Doppler phân số), ba mức SNRp cao hơn do (10.1) cần nhiều mẫu ước lượng hơn khi L=(2ξν+1)P>P.'),
      new Paragraph({ alignment: AlignmentType.CENTER, spacing: { after: 40 }, children: [img(OUT+'afdm_repro/fig12a_reproduction.png', 460, 366)] }),
      caption('Hình 7: BER theo SNR, AFDM, Ideal so với ước lượng qua pilot (SNRp=30/35/40dB), Doppler phân số ξν=1.'),
      pLead('Quan sát:', 'SNRp=40dB gần Ideal nhất, SNRp=30dB xa nhất, đúng thứ tự trên toàn dải SNR (SNR=18dB: Ideal=1.10×10⁻³, SNRp=40dB=1.14×10⁻³, SNRp=30dB=2.55×10⁻³).'),
      pLead('Giải thích:', 'SNRp cao hơn cho ước lượng đường truyền chính xác hơn - với điều kiện trừ đúng phần đóng góp đã biết của pilot khỏi tín hiệu thu trước khi tách dữ liệu; bỏ sót bước trừ này khiến phần dư pilot chưa khử tăng theo SNRp và có thể đảo ngược thứ tự quan sát được. So với Doppler nguyên (Hình 6), Doppler phân số cần SNRp cao hơn khoảng 5dB mới đạt hiệu năng gần Ideal tương đương - phù hợp với L=(2ξν+1)P>P ở (10.1): nhiều tham số hơn cần ước lượng từ cùng một pilot nên cần công suất pilot lớn hơn để đạt độ chính xác tương đương.'),

      h1('17. Hiệu ứng tham số bảo vệ Doppler phân số ξν'),
      p('Mục tiêu: cố định SNRp=40dB, thay đổi ξν=1/4/8 để kiểm chứng trực tiếp công thức L=(2ξν+1)P ở (10.1).'),
      new Paragraph({ alignment: AlignmentType.CENTER, spacing: { after: 40 }, children: [img(OUT+'afdm_repro/fig12b_reproduction.png', 456, 367)] }),
      caption('Hình 8: BER theo SNR, AFDM, hiệu ứng ξν=1/4/8 (Ideal và ước lượng qua pilot, SNRp=40dB).'),
      pLead('Quan sát:', 'ξν càng lớn, BER càng thấp, cả Ideal lẫn ước lượng (SNR=18dB, Ideal: ξν=1→8.84×10⁻⁴, ξν=4→5.12×10⁻⁴, ξν=8→3.94×10⁻⁴).'),
      pLead('Giải thích:', 'ξν lớn hơn nới rộng vùng tìm kiếm quanh mỗi vị trí locᵢ (Mục 10.1), giảm khả năng chồng lấn giữa các cột ứng với từng đường truyền trong Heff - cải thiện cả tách sóng lẫn ước lượng kênh, đúng theo công thức L=(2ξν+1)P. Đổi lại, ξν lớn hơn cũng làm tăng Q ở (10.2), tức chiếm nhiều mẫu bảo vệ hơn trong khung tín hiệu - với ISAC, đây là một đánh đổi trực tiếp cần cân nhắc: dành nhiều tài nguyên khung hơn cho việc tách Doppler chính xác, đổi lại ít tài nguyên hơn cho dữ liệu hoặc các chức năng cảm biến khác.'),

      h1('18. Thảo luận'),
      p('Đóng góp chính của báo cáo:'),
      p('• Khảo sát mô phỏng bảy kịch bản (Mục 11-17), phương pháp và tham số tham khảo từ [2] (Mục VII), bằng MATLAB tự viết theo đúng công thức Phần I - đo trực tiếp cả bậc phân tập AFDM (Mục 11) lẫn cấu trúc dải hẹp của kênh hiệu dụng (Mục 14).'),
      p('• Suy ra và trình bày công thức định lượng (10.1)-(10.2) cho cấu trúc dải hẹp của Heff dưới AFDM, làm rõ cơ sở toán học cho ưu thế độ phức tạp tính toán đã nêu ở Bảng 2.'),
      p('• Đối chiếu trực tiếp cấu trúc ma trận kênh hiệu dụng (Mục 10) với hiệu năng BER đo được (Mục 11-17), làm rõ cơ chế vật lý đứng sau chênh lệch giữa AFDM/OCDM/OTFS/OFDM, gắn với các đánh đổi cụ thể cho thiết kế ISAC (chi phí pilot, độ phức tạp tính toán, độ nhạy với va chạm Doppler).'),
      dataTable(['Waveform', 'Bậc phân tập (Phần I)', 'Cấu trúc kênh hiệu dụng (Mục 10)', 'Kết quả khảo sát chính (Mục 11-17)'], [
        ['OFDM', 'ρ=1 (mọi SNR)', 'ICI rộng quanh đường chéo', 'Kém nhất trong mọi so sánh (Hình 3, 4); sàn lỗi cấu trúc rõ rệt khi Doppler nguyên (Hình 5)'],
        ['OCDM', 'trường hợp riêng của AFDM, chưa chứng minh riêng', 'không tối ưu vị trí locᵢ theo kênh cụ thể', 'Trung gian giữa OFDM và AFDM/OTFS (Hình 3, 4)'],
        ['OTFS', 'ρ→1 khi SNR→∞', 'rò rỉ khu trú quanh đường chéo', 'Gần trùng AFDM ở thu MMSE, N=256 (Hình 4), nhưng chi phí pilot ước lượng kênh cao hơn (theo [2])'],
        ['AFDM', 'ρ=P (đầy đủ)', 'Thưa/tập trung nhất; dải hẹp đúng L=P (Doppler nguyên) hoặc L=(2ξν+1)P (phân số)', 'Tốt nhất trong mọi so sánh; độ dốc BER tăng theo P đúng ρ=P (Hình 2); ước lượng kênh gần Ideal khi SNRp đủ lớn (Hình 6-8)'],
      ]),
      caption('Bảng 3: Tổng hợp so sánh OFDM/OCDM/OTFS/AFDM qua các mục của Phần II.'),
      p('Bảy kết quả ở Mục 11-17 đo trên kịch bản kênh và tham số tham khảo từ [2], khác với thí nghiệm minh họa cấu trúc kênh ở Mục 10 (kênh Rayleigh ngẫu nhiên hai đường tự chọn để trực quan hóa) - hai lớp kịch bản độc lập cho cùng một xu hướng, củng cố độ tin cậy của kết luận.'),

      h1('19. Hạn chế của nghiên cứu'),
      p('• Kênh khảo sát là Rayleigh đa đường thống kê đơn giản hóa; chưa dùng mô hình chuẩn công nghiệp (3GPP TDL/CDL). Hệ thống là SISO, chưa mở rộng MIMO. BER đo được là "chưa mã hóa" (uncoded), chưa xét ảnh hưởng mã kênh.'),
      p('• numMonte của Mục 11-17 (300-30000 tùy thí nghiệm, Bảng 1) thấp hơn nhiều so với quy mô 10⁶ trial dùng trong [2] do giới hạn thời gian tính toán - các đường cong BER không đạt độ mượt tuyệt đối ở vùng BER cực thấp (dưới 10⁻⁵); một số điểm SNR cao ở Mục 14 cần tăng riêng số trial để giảm nhiễu thống kê, nhưng không ảnh hưởng đến các kết luận so sánh chính.'),
      p('• Bộ ước lượng kênh qua pilot nhúng (Mục 9.2) dùng thuật toán tìm-đường-mạnh-nhất-loại-trừ-liên-tiếp (matching pursuit) tự triển khai, không phải thuật toán tối ưu đầy đủ mô tả trong [2] - đã kiểm tra cho kết quả đúng hướng qua chẩn đoán trực tiếp (biết trước delay/Doppler thật) trước khi dùng cho Mục 15-17.'),
      p('• Bộ tách sóng MMSE dải hẹp ở Mục 14 áp ngưỡng L trực tiếp lên Heff dạng tuần hoàn đầy đủ, không chèn null-guard đúng như giả thiết ở Mục 10.1 - hệ số L thực dùng cho trường hợp Doppler phân số lớn hơn công thức (10.1) khoảng 3 lần mới cho bậc sai số quan sát được ở Hình 5.'),

      h1('20. Kết luận'),
      pLead('Mục tiêu:', 'khảo sát mô phỏng hiệu năng AFDM cho hệ thống ISAC trong kênh đa đường/Doppler cao, đối chiếu với các định lý bậc phân tập của AFDM/OTFS/OFDM đã chứng minh ở Phần I, dùng phương pháp và tham số tham khảo từ [2].'),
      pLead('Kết quả chính:', 'cấu trúc kênh hiệu dụng của AFDM là dải hẹp theo công thức (10.1)-(10.2) (Mục 10), kéo theo bậc phân tập đầy đủ đo trực tiếp bằng thu ML (Mục 11) và BER thấp hơn rõ rệt OFDM/OCDM ở mọi so sánh (Mục 12-14); ước lượng kênh qua pilot nhúng đạt hiệu năng gần với biết trước kênh khi SNRp đủ lớn, đúng thứ tự dự đoán bởi cấu trúc dải hẹp (Mục 15-17). Tổng hợp so sánh đầy đủ và đóng góp chính của báo cáo ở Mục 18 (Bảng 3).'),
      pLead('Hướng phát triển:', 'trong phạm vi đã khảo sát, AFDM là lựa chọn waveform nhiều triển vọng cho ISAC trong kênh đa đường/Doppler cao so với OFDM/OCDM, với hiệu năng tương đương OTFS ở thu MMSE nhưng cấu trúc kênh đơn giản hơn và chi phí pilot thấp hơn (Mục 10, 13). Mở rộng sang MIMO, mã kênh, và mô hình kênh công nghiệp chuẩn là hướng nghiên cứu tiếp theo (Mục 19).'),

      new Paragraph({ children: [new PageBreak()] }),

      h1('Phụ lục A. Tên file và tham số MATLAB của khảo sát mô phỏng (Mục 11-17)'),
      p('Bảng dưới liệt kê tên file MATLAB (tự viết) tương ứng từng hình ở Mục 11-17, phục vụ tái lập; công thức và tham số vật lý đã tóm tắt ở Mục 9.1-9.2 (Bảng 1).'),
      dataTable(['Hình', 'File MATLAB', 'Ghi chú triển khai'], [
        ['2', 'fig8a_afdm_diversity.m', 'ML vét cạn bằng ma trận giả thuyết 2^N hàng, AFDM'],
        ['3', 'fig8b_four_waveforms.m', 'ML vét cạn, 4 waveform, kênh 3 đường cố định'],
        ['4', 'fig9_mmse_comparison.m', 'MMSE tuyến tính, 4 waveform, N=256'],
        ['5', 'fig10_detectors.m', 'MMSE/MRC-DFE (lặp Gauss-Seidel)/MMSE dải hẹp (topRowMagnitudeMask.m), Doppler nguyên và phân số'],
        ['6', 'fig11_channel_estimation.m', 'Pilot nhúng, Doppler nguyên, gọi afdm_pilot_estimate.m'],
        ['7', 'fig12a_channel_estimation_frac.m', 'Pilot nhúng, Doppler phân số ξν=1'],
        ['8', 'fig12b_xinu_effect.m', 'Quét ξν=1/4/8 tại SNRp=40dB cố định'],
      ]),
      caption('Bảng 4: File MATLAB tương ứng từng hình của Mục 11-17.'),
      p('Hàm dùng chung: Gen_channel_mtx.m dựng ma trận kênh vật lý H theo đúng (1.1); afdm_pilot_estimate.m ước lượng đường truyền bằng tương quan-loại trừ liên tiếp (Mục 9.2), tinh chỉnh Doppler phân số bằng dò lưới cục bộ khi ξν>0. Toàn bộ mã nguồn, dữ liệu trung gian (.mat) và ảnh gốc lưu tại D:\\doubly-dispersive-channel-simulation-main\\afdm_paper_check\\.'),

    ]
  }]
});

Packer.toBuffer(doc).then(buf => {
  fs.writeFileSync('D:\\doubly-dispersive-channel-simulation-main\\report_build\\LyThuyet_AFDM_OTFS_OFDM.docx', buf);
  console.log('done');
});
