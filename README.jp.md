# fastapl -- 配列データの単発的な処理をPerl１行スクリプト風に行う。


## 始めに、
`% export PERL5LIB=$PERL5LIB:THIS_DIRECTORY/modules/perl5`
`% fastapl5 --help`

## 回帰テストを行う。
`% perl t/runRegressionTests.pl


## 要約
fastaplは断片的なPerlコードでfasta形式の配列データの様々な処理を便利にする。fastqplは、fastaplの姉妹プログラムでfastq形式のデータを扱う。fastaplはユーザが提供するコードを包んでPerlで実行するラッパですが、逆相補鎖の計算(rc)や塩基配列のアミノ酸配列への翻訳(translate, translations)など配列処理を便利にする関数も提供する。また、fastaplコマンド行で利用する地味なツールではあるが、xterm端末またはhtml形式で太字や色を使った表示はできる。例えば、color red, matches "CACA=1"でCACAに一致する配列領域を赤にして表示する。このmatchesのように配列領域を扱うfastaplの関数はDNAが二本鎖を持ち、場合によっては冠状であることを考慮して設計してあり、いろなケースに対応できる。

### 注意点
Linux環境下でv5.18.2などいくつかのPerlバージョンでの動作確認はしているが、それ以外の環境での検証は行われていない。


### 使用例

#### fastaplの使用例
* 長い配列を39残基目でちょんぎる。<BR>
`% fastapl5 -pe 'trim 0,38'  proteins.fasta`

* DNAの逆相補鎖を表示する。<BR>
`% fastapl5 -p  -e '$seq = reverse $seq; $seq =~ tr/acgtACGT/tgcaTGCA/'  shortDNA.fasta`

* Ｌで終る配列の割合を計算する。<BR>
`% fastapl5 -e '++$tot; ++$m if( $seq =~ /L$/ );'  -f 'say "$m/$tot"'  proteins.fasta`

* Ｍ以外の残基で始まる配列を表示する。<BR>
`% fastapl5 -ge '$seq !~ /^M/'  proteins.fasta`

* ふたつの入力ファイルにある配列群の配列長の合計を計算する。<BR>
`% fastapl5 -e '$tot += length $seq'  --end 'print "total seq. length = $tot\n"' proteins1.fasta proteins2.fasta`

* 配列群のアミノ酸組成を計算する。<BR>
`% fastapl5 -e '++$tot, $a{$_}+=100 for @seq'  --end 'printf "%s: %5.2f\n", $_, $a{$_}/$tot  for sort keys %a'	proteins.fasta`

* 各配列の最初の残基を飛ばして、配列群のアミノ酸組成を計算する。<BR>
`% fastapl5 -e '++$tot, $a{$_}+=100 for @seq[1..fin]'  --end 'printf "%s: %5.2f\n", $_, $a{$_}/$tot  for sort keys %a'  proteins.fasta`

* 各配列のアミノ酸組成を計算する。<BR>
`% fastapl5 -e '%a=(); $a{$_}+=100 for @seq;  print $id;  printf " $_:%4.2f", $a{$_}/len  for sort keys %a;  say ""'  proteins.fasta`

* 各レコードの配列をランダムにシャッフルする。<BR>
`% fastapl5 -b 'srand 7'  -pe '$seq = join "", shuffle @seq'  proteins.fasta`

* シーケンサーエラーの可能性を考慮しながら、ポリＡ鎖を切り捨てる。<BR>
`% fastapl5 -p  -e '$s = 0; $p = 0; $m = -9; for $i (reverse 0..$#seq){ $s += $seq[$i] eq "a" ? 1 : -2; if( $s > $m ){ $m = $s; $p = $i }} $m < 8 or substr( $seq, $p ) = ""'  DNAwithPolyAtail.fasta`

* 配列レコードをID順に並べ替える。<BR>
`% fastapl5 --sort -e '$id1 cmp $id2'  proteins.fasta`

* 配列を配列長に沿って並べ替える。<BR>
`% fastapl5 -se 'len1 <=> len2'  proteins.fasta`

* 最大100文字の行長で配列を表示する。<BR>
`% fastapl5 -pl 100  proteins.fasta`

* 最初の３つの配列レコードを表示する。<BR>
`% fastapl5 -ge '$c++ < 3'  proteins.fasta`

* 指定したIDの一致する配列レコードを表示する。<BR>
`% fastapl5 -b '%keep = qw(Q12154|GET3_YEAST 1 P24583|KPC1_YEAST 1 P38179|ALG3_YEAST 1)' -ge '$keep{$id}'  proteins.fasta`

* id.txtにあるIDに一致する配列レコードを表示する。(ids.txtの内容は１行にひとつのID)。<BR>
`% fastapl5 -b '%keep = map {chop;$_,1} <STDIN>' -ge '$keep{$id}'  proteins.fasta  <  ids.txt`

* (長い配列でも)１本の配列が１行に表示される形式で配列レコードを表示する。<BR>
`% fastapl5 -p -1  proteins.fasta`

* 配列レコードのhead行('>'から始まる行)でIDの次のフィールドに'mito'のあるレコードを表示する。<BR>
`% fastapl5 -ge '$head[1] =~ /mito/'  proteins.fasta`

* 各配列の順鎖と逆相補鎖をぞれぞれ表示する。逆相補鎖のIDに'_Crick'を付ける。<BR>
`% fastapl5 -pe '$id .= "_Crick"; pr; rc; $id =~ s/Crick$/Watson/'  shortDNA.fasta`

* DNA順鎖をRNAの転写配列に変換する。<BR>
`% fastapl5 -pe '$seq =~ tr/tT/uU/'  shortDNA.fasta`

* ６通りの考えられる読み枠をアミノ酸配列に仮想翻訳して表示する。<BR>
`% fastapl5 -e 'say $id;  say "  $_"  for translations'  DNAwithPolyAtail.fasta`

* 配列長が29アミノ酸以上のORF断片を配列レコードのコメントに書き込む。<BR>
`% fastapl5 -pe '$comment .= sprintf "# %2d %s\n", @$_[2,3]  for ORFs 29'  DNAwithPolyAtail.fasta`

* 配列長が29アミノ酸以上のORFを配列レコードのコメントに書き込む。<BR>
`% fastapl5 -pe '$comment .= sprintf "# %2d %s\n", @$_[2,3]  for ORFs "m29"'  DNAwithPolyAtail.fasta`

* 配列を大文字にして表示する。<BR>
`% fastapl5 -pe '$seq = uc;'  DNAwithPolyAtail.fasta`

* 各配列の２塩基頻度を表示する。<BR>
`% fastapl5 -e 'say hstr kmers 2'  mitoGenomes.fa`

* 各配列の２塩基頻度をコンパクトな形式で表示する。<BR>
`% fastapl5 -e 'say  "$id\t",  hstrf ": ", kmers 2'  mitoGenomes.fa`

* 配列群の２塩基頻度を表示する。<BR>
`% fastapl5 -e 'addKmers %c, 2' -f 'say hstr %c'  mitoGenomes.fa`

* 20塩基長以上の単塩基の繰り返し配列を表示する。<BR>
`% fastapl5 -e '$m = runs 20;  say  join "\t",  $a, $b, $b-$a, seq $a  while $m->ab'  yeast_chrII.fa`

* CACAに一致する個所を赤で表示する。逆相補鎖の一致個所には下線も引く。<BR>
-H -pe 'color red, matches "CACA=1"; ul matches "CACA_1"'  shortDNA.fasta

* 配列を１行100文字の形式で配列レコードを表示する。<BR>
`% fastapl5 -pl 100  proteins.fasta`

* 配列レコードの先頭行('>'から始まる行)でIDの次のフィールド'E.R.'のあるレコードを数える。<BR>
`% fastapl5 -ce '$head[1] eq "E.R."'  proteins.fasta`

* 配列レコードの先頭行('>'から始まる行)に記載してある情報をキー値ペアとして表示する。<BR>
`% fastapl5 -e 'say "  $id"; say hstr asHash'  proteinsHashHeaders.fasta`

* 芳香族を太字にし、目盛付きで配列を表示する。<BR>
`% fastapl5 -r -pe 'bold matches "[YWF]"'  proteins.fasta`

* DNAを二重鎖として表示する。<BR>
`% fastapl5 -P  shortDNA.fasta`

* 各配列の最後の４残基に色を付ける。<BR>
`% fastapl5 -pe 'colorize -4,0'  proteins.fasta`

* 各配列を目盛付きで二本鎖として表示する。配列に並べ、色付きで開読み枠を翻訳したトラックも表示する。<BR>
`% fastapl5 --ORF a110m80c  -Pr  -l100  InfluenzaA_H1N1_Iwate200903.fa`

* 各配列の(ゼロから数えて)5000塩基目から100塩基の領域を目盛付きで出力する。<BR>
`% fastapl5 -Re 'say "$id 5000,5100"; ps 5000,5100'  mitoGenomes.fa`

* 残基に色を付けて配列を表示する。<BR>
`% fastapl5 -C  proteins.fasta`

* DNA配列として、塩基に色を付けて配列を表示する。<BR>
`% fastapl5 -Cd  shortDNA.fasta`

* 配列レコードの先頭行('>'から始まる行)の最後のフィールドを出力する。<BR>
`% fastapl5 '-F;' -e 'say $head[-1]'  semicolonHeadFields.fasta`

* 配列を環状のDNA分子と見なし、原点を中心に10塩基長の逆相補鎖領域を表示する。<BR>
`% fastapl5 -Oe 'say seg 5, -5'  mitoGenomes.fa`

* 順鎖のTATAに一致する個所を赤にする。'%'は重複するマッチを含む意味(例えば...TATATA...)のTATATAは全部赤になる。<BR>
`% fastapl5 -pe 'color red, matches "TATAT%"'  m.fa`

* 二本鎖のどちらの鎖でもGAATATが一致する個所を白黒反転で表示する。<BR>
`% fastapl5 -OPe 'iv matches "GAATAT="'  mitoGenomes.fa`

* (配列を格納する)$seqの文字を小文字にしてから、tata[at]aの一致箇所を大文字にする。DNA配列を示す-dオプションでは'w'は[at]として扱われる。<BR>
`% fastapl5 -dpe '$seq = lc; upcase matches "tatawa"'  yeast_chrII.fa`

* アデニンの後に２アミノ酸「ＭＡ」に対応する２コドンが続く一致箇所を出力する。拡張正規表現を利用した例。<BR>
`% fastapl5 -db 'say expandRegex "(?i)a<MA>"'`

* 開始コドンを数えて出力する。入力配列はRNAでもDNAでも(大文字でも小文字)使える。<BR>
`% fastapl5 -e '$c=()=/(?i)(A[UT]G)/g; say "$id $c"'  yeast_chrII.fa`


## fastqplの使用例
* idが'4401B'を含むレコードを表示する。<BR>
`fastqpl5 -g  -e  '$id =~ /4401B/'  longreads_original_sanger.fastq`

* 各レコードに対し、IDとすべての塩基判定が正しい確立を出力する。<BR>
`fastqpl5 -e '$p = 1;  $p *= (1 - $_) for( errProb() ); print "$p\n"'  longreads_original_sanger.fastq`

* fastq形式をfasta形式に変換する。<BR>
`fastqpl5 -e 'print ">$head\n$seq\n"'  longreads_original_sanger.fastq`

* 平均Phred誤読確率スコアの降順にレコードを並べ替える。<BR>
`fastqpl5 -se '(sum Phred1) / length $seq1  <=>  (sum Phred2) / length $seq2'  longreads_original_sanger.fastq`

* 各配列の最初の４塩基を切り捨てる。<BR>
`fastqpl5 -pe 'trim 4'   longreads_original_sanger.fastq`

* 別な方法で、各配列の最初の４塩基を切り捨てる。<BR>
`fastqpl5 -pe 'trim 0,-5'   longreads_original_sanger.fastq`

* 配列とPhred誤読確率スコアを、'FSRRS4401BE7HA  t:37 c:37 a:37 g:35 T:35 T:35...'のような形式で出力する。<BR>
`fastqpl5 -e '@p = Phred; @v = pairwise {"$a:$b"} @seq, @p; say "$id  @v"'   longreads_original_sanger.fastq`

* 各配列の後端にある、誤読確率が0.1以上の塩基を切り捨てる。<BR>
`fastqpl5 -ge  'trim  0,  lastidx {$_ < 0.1} errProb;  $seq'   longreads_original_sanger.fastq`

* 各配列に最大30塩基長の擬似ポリＡ鎖を付け加える。ポリＡ鎖の誤読確率が塩基毎に増えるようにする。<BR>
`fastqpl5 -b 'srand 5'  -pe  '$len = rand 30; $p = .05; @p = map {$p *= 1.1} 1..$len; $seq .= "a" x $len; $quality .= prob_to_qualStr @p'  longreads_original_sanger.fastq`

* BWA trimmingのような処理を行う。調整クオリティスコア(誤読確率Phredスコアから15引いた値)の和が最大の接頭配列(同点の場合は短い接頭配列が優先)を出力する。<BR>
`fastqpl5 -ge '$sum = 0;  ($i,$max) = maxidx apply {$_ = $sum += $_-15} Phred;  trim 0,$i;  $max > 0'  longreads_with_dummy.fastq`

* Solexa trimmingのような処理を行う。調整クオリティスコア(誤読確率Phredスコアから15引いた値)の和が最大の部分配列を出力する。<BR>
`fastqpl5 -ge '($beg,$end,$max) = maxSeg apply {$_ -= 15} Phred;  trim $beg,$end;  $max'  longreads_with_dummy.fastq`


### 関連ツール
perltab
