# GROMACSチートシート

## 計算準備編

### ndxファイルを作る

`gmx make_ndx -f initial.gro -o index.ndx`

この時コマンド内でもグループを作る方法が色々ある。
例を載せておく。
「a P」とすることで、全ての"P"という原子を選択できる。
「1 | 2 | 3」とすることで、グループ1または2または3に含まれている原子が選択できる。
これにより、これらがi番目のグループとして登録される。
また名前を変えるには、「name i(グループ番号) Phosphorous(つけたいグループ名)」とすればよい
あとは"q"を押してセーブする。

またEnterを押せば登録されたグループのリストを確認できる。

さらに-nオプションで既存のindexファイルを読み込んで新たなindexファイルを作ることもできる。

### 系を複製したいとき(groファイルだけでなくpdbファイルも可)

`gmx genconf -f input.gro -o double.gro -nbox 2 1 1`

これで系をx軸方向に2倍、y,z軸方向に1倍ずつ、つまり元の系を2倍に複製できる。
pdbでもできるらしいので便利。
これ周期境界設定していないと同じ位置に重なってうまく複製できないので注意！！！！！！！！

### groファイルの残基番号を修正したいとき

`gmx editconf -f initial.gro -o number.gro -resnr 1`

-resnrは何番から付け直すかを指定できる。

但しうまく直らない場合もあるのでその都度生成したファイルの中身を確認すること。



## 計算編

### 速度テストしたいとき

`gmx mdrun -deffnm npt -nsteps 1000`

とすると、強制的に1000stepで終了できる。
いちいち.mdpファイルを編集しなくていいから助かる。

### 途中で終了してしまったrunの再開

`gmx mdrun -deffnm npt -cpi npt.cpt`

npt.cptはチェックポイントファイルでありMD計算の一定間隔で状態を記録・updateしている。
そのためこれを使ってリスタートすると、このチェックポイントファイルの状態からスタートする。
ちなみに`gmx check -f npt.cpt`でどのステップのものか確認できるらしい！
便利すぎる...
ちなみにnpt_prev.cptは最新の1つ前のcptだと思われる。



## 解析準備編

### 周期境界条件を直す

`gmx trjconv -s npt.tpr -f npt.xtc -o npt_noPBC.xtc -pbc mol`

必要に応じて-centerオプションを使う。

### トラジェクトリのフレームを落とす

トラジェクトリのフレームを落としたい時は-skipオプションで以下のようにする。

`gmx trjconv -f npt_noPBC.xtc -o npt_cut.xtc -skip 10`

これで10フレーム間隔でトラジェクトリを少なくできる。

### トラジェクトリを結合したいとき

`gmx trjcat -f npt1.xtc npt2.xtc npt3.xtc -o npt_all.xtc -settime`

これでnpt1.xtcからnpt3.xtcまでを全て結合できる。
この時settimeとしてひたすらにcを入力するあるいはcatとしないと、トラジェクトリを連結してくれない。

なんかトラジェクトリを重ねてしまって(?)、生成後のトラジェクトリのフレーム数も一致しない。

この辺りは少しよく分かっていない。

## 解析編

### 系の物理量を調べる

`gmx energy -f npt.edr -o pressure.xvg`

例えばこのコマンドを実行した後に18 0(バージョンによってナンバー違う)を入力すると時間に対しての圧力が出力できる。
圧力のみならず、温度、密度、ポテンシャルエネルギー、系の全ての合計のエネルギー(total energy)、さらには熱浴も込みでの合計のエネルギー(conserved energy)、ボックスのz軸方向の長さなども調べれてしまう。
このコマンドはやばい。
あと拘束のエネルギー(restraint potential)も調べれるはずだけど、特に距離に関してはbondsで拘束していたのでこのコマンドで調べれなかった。
distant restraintを使っていたらこれで調べれていたかもしれない。

### 拡散係数を調べる

まず、計算前に周期境界で変になる結合を直したほうがいい。

実際のコマンドとしては
`gmx msd -f npt.xtc -s npt.tpr -o msd.xvg (-tu ns)`

ちなみに例えば500 psの計算であればその内の10%と90%すなわち50ps-450psがmsd計算に使用される。
これはそれぞれ-beginfitと-endfitで使用したい時刻を決めれる。
またmsdの計算ではどれくらいの時間幅を考えるかを決める必要があり、デフォルトでは10 psとなっている。
理想を言えば500 psを1 ps刻みでやりたいが、計算コストがかかるのとそこまで時間幅を短くする意味もない。
そのためこの場合は0,10,20,...,490,500という時間幅の平均二乗変位をプロットすることになる。

### 動径分布を調べる

`gmx rdf -f npt.xtc -s npt.tpr -n rdf.ndx -o rdf.xvg`

ちなみに-exclオプションを使うと、結合している原子同士のrdfを考えなくなるっぽい。つまりOH2をrefにしてH1とH2をselにすると
何のピークもない結果が得られるらしい。
計算時間が少なくて短いからか、ネットのやつとは一致しなかったけど多分使用の仕方は問題なさそう。

但しgroファイルを使ってrdfを求めようとすると、原子の質量や原子半径が分からないために残基や原子名から推測して使用されるらしい。
これはlammpsで粗視化計算したものを無理やりgroにした場合にはこの方法は使えないということ。またトラジェクトリに結合が変になってるやつがあるとまずいって。そうなるとやっぱりlammpsはlammpsでrdf求める必要がある。

追記: 質量に関してはpSPICAでもpsfファイルに書いてあるから問題ない。但しlammpsだと不便すぎるのでMPDynを使うことにした。
また細かいオプションについても理解したので記しておく。

`gmx rdf -f npt.xtc -s npt.tpr -o rdf.xvg -n rdf.ndx -rmax 8 -bin 0.1`

-rmaxはrdfをどの程度の距離まで測るか。この場合は8 nmまで測り、またその測る刻みを0.1 nmにするということ。

### 二次元密度マップを調べる

動径分布は距離の情報を与えてくれはするが、分子がどの位置に集まるとかは当然のことながら分からない。

そういう場合は二次元(x,y平面)でマップを作れば一目瞭然である。コマンドは以下。

`gmx densmap -f npt.xtc -s npt.tpr -n analy.ndx -o map.xpm -dmax 2`

xpmという形式でマップを出力する。-dmaxは密度(グレースケール)の最大値。

さらにこのxpmを**xpm2ps**コマンドを使うことでepsに変換しつつ、マップを綺麗に整えることが可能。コマンドは以下。

`gmx xpm2ps -f map.xpm -di input.m2p -o map.eps -title ylabel -rainbow blue`

blueにすると密度の低い方から高い方を青から赤にできる。逆にしたければ-rainbow redにする。

さらにこのコマンド実行時には**input.m2p**というインプットを読み込む。これはグラフの設定をするファイルで、gnuplotのような

イメージをもっておくと分かりやすい。インプットの中身は次のように設定している。

```
; Command line options of xpm2ps override the parameters in this file
black&white              = no           ; Obsolete
titlefont                = Times-Roman  ; A PostScript Font
titlefontsize            = 10           ; Font size (pt)
legend                   = yes          ; Show the legend
legendfont               = Times-Roman  ; A PostScript Font
legendlabel              =              ; Used when there is none in the .xpm
legend2label             =              ; Used when merging two xpm's 
legendfontsize           = 30.0           ; Font size (pt)
xbox                     = 1         ; x-size of a matrix element
ybox                     = 1          ; y-size of a matrix element
matrixspacing            = 20.0         ; Space between 2 matrices
xoffset                  = 100.0          ; Between matrix and bounding box
yoffset                  = 100.0          ; Between matrix and bounding box
x-major                  = 20           ; Major ticks on x axis every .. frames
x-minor                  = 5            ; Id. Minor ticks
x-firstmajor             = 0            ; First frame for major tick
x-majorat0               = no           ; Major tick at first frame
x-majorticklen           = 0          ; x-majorticklength
x-minorticklen           = 0          ; x-minorticklength
x-label                  = X             ; Used when there is none in the .xpm
x-fontsize               = 0.1           ; Font size (pt)
x-font                   = Times-Roman  ; A PostScript Font 
x-tickfontsize           = 0.1          ; Font size (pt)
x-tickfont               = Helvetica    ; A PostScript Font
y-major                  = 20
y-minor                  = 5
y-firstmajor             = 0
y-majorat0               = no
y-majorticklen           = 0
y-minorticklen           = 0
y-label                  = Y
y-fontsize               = 0.1
y-font                   = Times-Roman
y-tickfontsize           = 0.1
y-tickfont               = Helvetica
```

但しこの設定はけっこう融通が効かない。例えばラベルとグラフの位置を調節できないので、文字の大きさを大きくしていくとグラフと被る。そのためラベルの大きさを0.1にすることでほとんどラベルを見えなくして、パワポなどに載せる時にテキストボックスで自分で書いている。またメモリを付けたければtickfontsizeをいじればいいが、つけるとシンプルに汚くなる。そのためこれも同様に0.1にして見えなくする。ちなみにメモリを消すというオプションはないようである。また0はプログラムに認識されないため0.1にする。

### 傾斜角を調べる

タンパク質のN末端からC末端へ向かうベクトルとz軸とのなす角を調べたかったので調べていたところこのコマンドを発見した。

`gmx gangle -f npt.xtc -s npt.tpr -n analy.ndx -oall tiltangle.xvg -g1 vector -g2 z -tu ns`

-g1 でベクトルと定義し、コマンド実行後にどのグループを用いるかを選ぶ。例えばanaly.ndxの中に[ PROTEIN1 ]というグループが登録されており、このリストでは1番(N末)と54番(C末)のみが定義されている。このグループを選ぶと1番の粒子から54番の粒子へ向かうベクトルが定義される。また-g2ではz軸オプションを指定しているためz軸方向のベクトルを定義できる。そしてこの2つのベクトルのなす角を計算することができる。-oavにすると時間経過の累積平均になるので、-oallにすることで実際のトラジェクトリ内での角度を出力できる。

###  ４つの原子からなる二面角を調べる

`gmx angle -f step6.1_equilibration.xtc -n dihedral.ndx -all -type dihedral -ov dihedraltest.xvg`

-f は入力ファイル、-nは4つの原子を指定するために必要なインデックスファイル、-allは(そのグループの？)平均値と別に個々の二面角を出してくれるらしいけど、今は４つの原子しかないから(４つの原子を１グループに.ndxでしている)同じ値が２列表示されてる。
-allオプション外すと一列しか表示されなかった。
-typeでどういう角度を調べるかを変えれる！これがないとデフォルトで結合角になるらしく、４つの原子をインデックスで指定してても３つの原子がないとか言われる。
なので-type dihedralとする必要がある。
もし指定しなければ結合角を調べることができるはず。
-ovは出力に関するオプションで、何を求めるかによって変わるが、単なる角度の出力なら-ovでいいはず

### ２つの原子の距離を調べる

`gmx distance -f step6.6_equilibration.xtc -n distance.ndx -oav distave.xvg`　(ndxで見たいものを別々のグループに記述する必要があるらしい)でいけるはず。

このコマンド結局使ったかどうか忘れたからちょっと怪しいかもしれない。
