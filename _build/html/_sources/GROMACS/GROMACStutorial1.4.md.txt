# 4.【補足】トポロジーファイルの詳細

ここでは`gmx pdb2gmx`で生成したトポロジーファイル（topol.top）について、その詳細を確認する。 構造ファイル（1AKI_processed.gro）には系の全原子の位置（座標）が記述されているのに対し、トポロジーファイルには原子の種類や結合の情報などが含まれている。 ただし、トポロジーファイルは系に分子を追加する際に更新されていくため、ここでは最終的なファイルを見ていくものとする。 （はじめに生成したものは水やイオンの情報が不足しているはず）

## 力場パラメータ、分子タイプセクション

まず、冒頭のいくつかのコメント行の後、力場パラメータに関する記述がある。 今回は力場にOPLS-AAを選択したので、これ以降のパラメータはOPLS-AAで導入されるものとなっている。

```
topol.top
; Include forcefield parameters
#include "oplsaa.ff/forcefield.itp"
```

その次は[ moleculetype ]となっていて、タンパク質の名前と相互作用パラメータの扱い方が記述されている。 ここでは原子間の相互作用に関して、結合が3原子以内のものは非結合性相互作用を無視するという設定をしている。 これ以上の説明はここではしないので、詳細はGROMACSのマニュアルを読んでみてほしい。

```
topol.top
[ moleculetype ]
; Name            nrexcl
Protein_chain_A     3
```

## 原子セクション、結合関連セクション

続いての[ atoms ]にはタンパク質の構成原子の情報が書かれている。 各原子に1から番号が割り振られ、それがタンパク質で何番目のどの残基に属する原子であるか、元素の種類、電荷、質量などが記述されている。 次の例では、タンパク質の1番の原子は1番目のLys残基に属する窒素原子であり、電荷-0.3、質量14.0027であることがわかる。

```
topol.top
[ atoms ]
;   nr       type  resnr residue  atom   cgnr     charge       mass  typeB    chargeB      massB
; residue   1 LYS rtp LYSH q +2.0
     1   opls_287      1    LYS      N      1       -0.3    14.0027   ; qtot -0.3
     2   opls_290      1    LYS     H1      1       0.33      1.008   ; qtot 0.03
     3   opls_290      1    LYS     H2      1       0.33      1.008   ; qtot 0.36
     4   opls_290      1    LYS     H3      1       0.33      1.008   ; qtot 0.69
     5  opls_293B      1    LYS     CA      1       0.25     12.011   ; qtot 0.94
     6   opls_140      1    LYS     HA      1       0.06      1.008   ; qtot 1
```

また、続く[ bonds ]、[ angles ]、[ dihedrals ]では、タンパク質の結合、結合角、二面角に関する情報が記述されている。 ここでどの原子とどの原子の間に結合があるか、結合した結果どの位置に結合角や二面角が生じるかを知ることができる。

```
topol.top
[ bonds ]
;  ai    aj funct       c0       c1       c2       c3
    1     2     1 
    1     3     1 
    1     4     1 
    1     5     1 
    5     6     1 
    5     7     1 

[ pairs ]
;  ai    aj funct       c0       c1       c2       c3
    1     8     1 
    1     9     1 
    1    10     1 
    1    24     1 
    1    25     1 
    2     6     1 

[ angles ]
;  ai    aj    ak funct       c0       c1       c2       c3
    2     1     3     1 
    2     1     4     1 
    2     1     5     1 
    3     1     4     1 
    3     1     5     1 
    4     1     5     1

[ dihedrals ]
;  ai    aj    ak    al funct       c0       c1       c2       c3       c4       c5
    2     1     5     6     3 
    2     1     5     7     3 
    2     1     5    23     3 
    3     1     5     6     3 
    3     1     5     7     3 
    3     1     5    23     3 
```

なお[ pairs ]では分子内のペア相互作用を定義し、1-4非結合相互作用などを扱うことができる。 また[ dihedrals ]ではimproper dihedralも定義されているため、セクションが2つに分けられている。

## 位置拘束パラメータ、水・イオンの力場パラメータ

このセクションでは位置拘束パラメータの定義と、水とイオンの力場パラメータの導入をしている。 タンパク質の位置拘束パラメータは`gmx pdb2gmx`を実行した時に生成したposre.itpの中に記述されており、位置拘束を有効にする設定をした場合に呼び出される。 同様に、水の位置拘束も設定を有効にするとPOSRES_WATERから呼び出される。

また、水とイオンのパラメータもOPLS-AAで導入するので、SPC/E水モデルのパラメータファイルspce.itpとイオンのパラメータファイルions.itpを読み込んでいる。

```
topol.top
; Include Position restraint file
#ifdef POSRES
#include "posre.itp"
#endif
    
; Include water topology
#include "oplsaa.ff/spce.itp"
    
#ifdef POSRES_WATER
; Position restraint for each water oxygen
[ position_restraints ]
;  i funct       fcx        fcy        fcz
   1    1       1000       1000       1000
#endif
    
; Include topology for ions
#include "oplsaa.ff/ions.itp"
```

今回は水モデルにSPC/Eモデルを使用したが、他にもSPCやTIP3P、TIP4Pなどが選択できる。 水モデルにはそれぞれの特徴があるので、シミュレーションの目的に合ったモデルを選択するようにしよう。

## 系全体の情報

最後に、系の名前と含まれる分子の一覧が記述されているセクションがある。 この系での系の名前は「LYSOZYME」となっているが、この部分は任意の名前に変更することができる。 一方で、[ molecules ]には系に含まれる分子の名前と個数が記されているので、作成した系が正しいかどうか確認しておく。

```
topol.top
[ system ]
; Name
LYSOZYME
    
[ molecules ]
; Compound      #mols
Protein_A         1
SOL           10636
CL                8
```

以上がトポロジーファイルの内容である。 基本的には`gmx pdb2gmx`で生成されるものを用いればよいが、シミュレーションによっては一部修正をすることもあるだろう。
