// Unpack the BAM-encoded sequence data, which is of course not
// necessarily NUL-terminated.
string alignment::unpack_seq(const char* raw_seq, int seq_length) {
  static const char decode[] =
    "===A=C=M=G=R=S=V=T=W=Y=H=K=D=B=N"
    "A=AAACAMAGARASAVATAWAYAHAKADABAN"  //  1 = A
    "C=CACCCMCGCRCSCVCTCWCYCHCKCDCBCN"  //  2 = C
    "M=MAMCMMMGMRMSMVMTMWMYMHMKMDMBMN"
    "G=GAGCGMGGGRGSGVGTGWGYGHGKGDGBGN"  //  4 = G
    "R=RARCRMRGRRRSRVRTRWRYRHRKRDRBRN"
    "S=SASCSMSGSRSSSVSTSWSYSHSKSDSBSN"
    "V=VAVCVMVGVRVSVVVTVWVYVHVKVDVBVN"
    "T=TATCTMTGTRTSTVTTTWTYTHTKTDTBTN"  //  8 = T
    "W=WAWCWMWGWRWSWVWTWWWYWHWKWDWBWN"
    "Y=YAYCYMYGYRYSYVYTYWYYYHYKYDYBYN"  // ...and similarly for columns w.r.t.
    "H=HAHCHMHGHRHSHVHTHWHYHHHKHDHBHN"  // the second character of each pair.
    "K=KAKCKMKGKRKSKVKTKWKYKHKKKDKBKN"
    "D=DADCDMDGDRDSDVDTDWDYDHDKDDDBDN"
    "B=BABCBMBGBRBSBVBTBWBYBHBKBDBBBN"
    "N=NANCNMNGNRNSNVNTNWNYNHNKNDNBNN"; // 15 = N

  const unsigned char* packed = reinterpret_cast<const unsigned char*>(raw_seq);

  int odd = seq_length % 2;
  string seq(seq_length + odd, 'X');

  string::iterator it = seq.begin();
  for (int i = 0; i < seq_length; i += 2) {
    int ndx = 2 * *packed++;
    *it++ = decode[ndx];
    *it++ = decode[ndx + 1];
  }

  if (odd)
    seq.erase(seq_length);

  return seq;
}
