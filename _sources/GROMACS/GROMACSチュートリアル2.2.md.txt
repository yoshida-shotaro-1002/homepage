# 2. 力場パラメータとトポロジーの準備

前回はアルコール分子の構造を作成し、それを使ってアルコール二重膜系を作成した。 しかし、この構造を用いてシミュレーションを行うにあたって、GROMACSでこの分子を認識させ、力場を与えてやる必要がある。 今回はGROMACSでの力場の扱い方を学ぶためにCHARMM力場を導入し、アルコール分子のトポロジーをGROMACSに認識させて、実際にシミュレーションが開始できるまでの準備をする。

1. CHARMM力場のダウンロード
2. CHARMMトポロジーの追加
3. 構造ファイル（groファイル）の作成

また、このページでもCHARMM力場のダウンロードのためにインターネット環境が必要である。 通信環境を確認しておこう。

## CHARMM力場のダウンロード

まず、シミュレーションに使用する力場を決定する。 チュートリアル①ではGROMACSにデフォルトでインストールされているOPLS-AA力場を使用したが、今回は新たな力場をインストールして使えるように設定する。

使用する力場はCHARMM36力場である。 [CHARMM Force Field Files](http://mackerell.umaryland.edu/charmm_ff.shtml)にアクセスし、**CHARMM36 Files for GROMACS**のセクションから最新版のCHARMM36力場ファイル（charmm36-xxx.ff.tgz）をダウンロードする。

<img src="https://shalad2.github.io/JellyPage/gromacs/images/charmm_ff.png" width=80%>

ダウンロードできたら、このファイルをGROMACSが認識できる場所で展開する。 GROMACSをインストールした場所によって名前は異なるかもしれないが、
`~/App/gromacs-XXX/gmx_install/share/gromacs/top`
にcharmm36-xxx.ff.tgzを移動し、ここで中身を展開するだけでCHARMM36力場が使用可能となる。

```
# ダウンロードしたファイルを移動
mv charmm36-xxx.ff.tgz ~/App/gromacs-XXX/gmx_install/share/gromacs/top/
cd ~/App/gromacs-XXX/gmx_install/share/gromacs/top/
# ファイルを展開
tar xvf charmm36-xxx.ff.tgz
```

`App`や`gmx_install`などは[GROMACSのインストール](https://shalad2.github.io/JellyPage/gromacs/gromacs_install.html)で紹介している場所なので、もし環境が異なるようであれば参考にしてほしい。

### CHARMMトポロジーの追加

CHARMM-GUIで作成したアルコール分子はオリジナルの分子であるため、GROMACSはデフォルトでこの分子を認識することができない。 そこで、GROMACSが認識できる分子一覧のファイルに今回作成したアルコール分子のトポロジーを与えてやることによって、この問題を解決できる。

アルコール分子のトポロジーファイルには原子情報と結合情報が含まれており、これは前回CHARMM-GUIで分子を作成した際に同時に生成している。 アルコール分子にLALという名前をつけたので、トポロジーファイルはCHARMM-GUIファイルの`lal`内の`lal.rtf`という名前になっているはずである。 テキストエディタで内容を確認すると、次のようになっている。

```
lal.rtf
# 原子情報
RESI lal            0.000 ! param penalty=   0.000 ; charge penalty=   0.000
GROUP            ! CHARGE   CH_PENALTY
ATOM C1     CG331  -0.270 !    0.000
ATOM C2     CG321  -0.177 !    0.000
ATOM C3     CG321  -0.183 !    0.000
ATOM C4     CG321  -0.180 !    0.000
ATOM C5     CG321  -0.180 !    0.000

# 結合情報
BOND C1   C2  
BOND C2   C3  
BOND C3   C4  
BOND C4   C5  
BOND C5   C6
```

しかし、このトポロジーファイルはGROMACSの形式ではないため（ここが一番面倒だが）形式を変換してやる必要がある。 少し長いがサンプルを以下に示すので、コピーして使ってほしい。

```
[ LAL ] ; original molecule made by charmm-gui
 [ atoms ]
        C1     CG331  -0.270   1
        C2     CG321  -0.177   2
        C3     CG321  -0.183   3
        C4     CG321  -0.180   4
        C5     CG321  -0.180   5
        C6     CG321  -0.180   6
        C7     CG321  -0.180   7
        C8     CG321  -0.180   8
        C9     CG321  -0.180   9
        C10    CG321  -0.180  10
        C11    CG321  -0.180  11
        C12    CG321  -0.180  12
        C13    CG321  -0.180  13
        C14    CG321  -0.180  14
        C15    CG321  -0.182  15
        C16    CG321   0.054  16
        O      OG311  -0.651  17
        H1     HGA3    0.090  18
        H2     HGA3    0.090  19
        H3     HGA3    0.090  20
        H4     HGA2    0.090  21
        H5     HGA2    0.090  22
        H6     HGA2    0.090  23
        H7     HGA2    0.090  24
        H8     HGA2    0.090  25
        H9     HGA2    0.090  26
        H10    HGA2    0.090  27
        H11    HGA2    0.090  28
        H12    HGA2    0.090  29
        H13    HGA2    0.090  30
        H14    HGA2    0.090  31
        H15    HGA2    0.090  32
        H16    HGA2    0.090  33
        H17    HGA2    0.090  34
        H18    HGA2    0.090  35
        H19    HGA2    0.090  36
        H20    HGA2    0.090  37
        H21    HGA2    0.090  38
        H22    HGA2    0.090  39
        H23    HGA2    0.090  40
        H24    HGA2    0.090  41
        H25    HGA2    0.090  42
        H26    HGA2    0.090  43
        H27    HGA2    0.090  44
        H28    HGA2    0.090  45
        H29    HGA2    0.090  46
        H30    HGA2    0.090  47
        H31    HGA2    0.090  48
        H32    HGA2    0.090  49
        H33    HGA2    0.090  50
        H34    HGP1    0.419  51
 [ bonds ]
        C1   C2
        C2   C3
        C3   C4
        C4   C5
        C5   C6
        C6   C7
        C7   C8
        C8   C9
        C9   C10
        C10  C11
        C11  C12
        C12  C13
        C13  C14
        C14  C15
        C15  C16
        C16  O
        C1   H1
        C1   H2
        C1   H3
        C2   H4
        C2   H5
        C3   H6
        C3   H7
        C4   H8
        C4   H9
        C5   H10
        C5   H11
        C6   H12
        C6   H13
        C7   H14
        C7   H15
        C8   H16
        C8   H17
        C9   H18
        C9   H19
        C10  H20
        C10  H21
        C11  H22
        C11  H23
        C12  H24
        C12  H25
        C13  H26
        C13  H27
        C14  H28
        C14  H29
        C15  H30
        C15  H31
        C16  H32
        C16  H33
        O    H34
```

このトポロジーを先程ダウンロード、展開したcharmm36-xxx.ffの中の**merged.rtp**というファイルの末尾に追加する。 merged.rtpにはGROMACSが認識できる分子の一覧が書かれているので、ここにトポロジーを追加することで新しい分子を扱うことが可能となる。

なお、既存のファイルを変更する時は必ずバックアップを取っておくようにしよう。

## 構造ファイル（groファイル）の作成

最後に、GROMACSが新しい力場と分子を正しく認識し、構造ファイルを作成できることを確認しておく。 以前のチュートリアルでも扱った通り、groファイルの作成には`gmx pdb2gmx`を使用する。

```
gmx pdb2gmx -f system.pdb -o system.gro
```

コマンドを実行すると、力場の選択画面になる。 これまでの設定が正しくできているとCHARMM36が選択できるはずである。

```
Select the Force Field:
From '/home/user/App/gromacs-XXX/gmx_install/share/gromacs/top':
1: AMBER03 protein, nucleic AMBER94 (Duan et al., J. Comp. Chem. 24, 1999-2012, 2003)
2: AMBER94 force field (Cornell et al., JACS 117, 5179-5197, 1995)
3: AMBER96 protein, nucleic AMBER94 (Kollman et al., Acc. Chem. Res. 29, 461-469, 1996)
4: AMBER99 protein, nucleic AMBER94 (Wang et al., J. Comp. Chem. 21, 1049-1074, 2000)
5: AMBER99SB protein, nucleic AMBER94 (Hornak et al., Proteins 65, 712-725, 2006)
6: AMBER99SB-ILDN protein, nucleic AMBER94 (Lindorff-Larsen et al., Proteins 78, 1950-58, 2010)
7: AMBERGS force field (Garcia & Sanbonmatsu, PNAS 99, 2782-2787, 2002)
8: CHARMM27 all-atom force field (CHARM22 plus CMAP for proteins)
9: CHARMM36 all-atom force field (July 2020)
10: GROMOS96 43a1 force field
11: GROMOS96 43a2 force field (improved alkane dihedrals)
12: GROMOS96 45a3 force field (Schuler JCC 2001 22 1205)
13: GROMOS96 53a5 force field (JCC 2004 vol 25 pag 1656)
14: GROMOS96 53a6 force field (JCC 2004 vol 25 pag 1656)
15: GROMOS96 54a7 force field (Eur. Biophys. J. (2011), 40,, 843-856, DOI: 10.1007/s00249-011-0700-9)
16: OPLS-AA/L all-atom force field (2001 aminoacid dihedrals)
```

`9`を選択して次へ進むと、使用する水モデルの種類を聞かれる。 前回はSPE/Cモデルを使用したが、今回はTIP3Pモデルを使用する。

```
Select the Water Model:
1: TIP3P	TIP 3-point, recommended, by default uses CHARMM TIP3 with LJ on H
2: TIP4P	TIP 4-point
3: TIP5P	TIP 5-point
4: SPC		simple point charge
5: SPC/E	extended simple point charge
6: None
```

`1`を選択し、何もエラーが出なければgroファイルが生成されているであろう。 また同時に、大量のトポロジーファイルと位置拘束ファイルが生成されているであろう。 これでは扱うファイルが増えてしまい、見た目もあまりよくないので、ここからはファイルを1つにまとめてみたいと思う。



現在`gmx pdb2gmx`をしたことでトポロジーファイルが自動生成されたが、PDBファイルに含まれているアルコール分子256個が別々のファイルに分けられてしまっている。 しかしながら、256個のアルコール分子は全て同じものであるので、わざわざファイルを分ける必要はなく、むしろ全て同じトポロジーファイルを利用することが望ましい。 従って、topol_Other.itp及びposre_Other.itpのみを残して、topol_OtherXXX.itpやposre_OtherXXX.itpなどは全て削除してしまおう。 （ただし、topol.topは残しておくこと）

```
# 生成したファイル（一部のみ抜粋）
posre_Other.itp       topol_Other.itp
posre_Other2.itp      topol_Other2.itp
posre_Other3.itp      topol_Other3.itp

posre_Other256.itp    topol_Other256.itp
topol.top
# posre_Other.itp, topol_Other.itp, topol.topのみを残して、他のitpファイルを全て削除
```

さらに、系全体のトポロジーファイルであるtopol.topを編集する。 「Include chain topologies」の項目に削除したファイル一覧があるので、topol_Other.itpのみを残して他を削除する。

```
topol.top (modified)
; Include chain topologies
#include "topol_Other.itp"
```

また、ファイル末尾の[ molecules ]のセクションでも同様の修正をして、Otherの分子数を256に変更する。

```
topol.top (modified)
[ moleules ]
; Compound      #mols
Other             256
SOL              7680
```

変更を保存して、次の5つのファイルが手元にあれば、ここまでの操作は完了となる。 なお、ここからsystem.pdbとposre_Other.itpは使用しないが、削除せずに残しておこう。

1. system.pdb
2. system.gro
3. topol.top
4. topol_Other.itp
5. posre_Other.itp

## まとめ

ここまで、アルコール二重膜の力場とトポロジーの設定を行い、GROMACSでシミュレーションを実行できるところまで辿り着いた。 今回の設定をまとめると次の通りである。

1. `CHARMM36力場`をダウンロードし、GROMACSで読み込めるようにした。
2. アルコール分子のトポロジーを`merged.rtp`に追加し、GROMACSで分子を扱えるようにした。
3. groファイルを作成し、系全体のトポロジーファイルを修正した。

次回は作成したgroファイルとトポロジーファイルを用いてシミュレーションを実行する。 シミュレーションの設定ファイルであるmdpファイルの詳細を学んでいこう。

[PREVIOUS](https://shalad2.github.io/JellyPage/gromacs/gromacs_tutorial2_1.html)
