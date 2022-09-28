# VMDチートシート

Tclの説明はありません。



## 画像の変換

`mogrify -fuzz 10% -density 2000 -format png *tga`

ImageMagickのmogrifyコマンドでtgaをpngにする。



## 構造作成

### 分子のオブジェクト生成

`set val [atomselect top "all"]`

全ての分子を選択している。

### 分子を平行移動

`$val moveby {0 0 10}`

10Å分全ての分子をz軸方向に移動。

### 分子を回転移動

`$val move [transaxis z 30]`

z軸に対して30°全ての分子を回転移動。

### 分子のセンタリング

```
proc center { selText } {
	 set total [atomselect top all]
     set ref [atomselect top $selText]
     set nf [molinfo top get numframes]
     puts [format "%i frames\n" $nf]
     # Center the reference system around (0, 0, 0)
     for {set f 0} {$f < $nf} {incr f} {
     $ref frame $f
     $total frame $f
     $total moveby [vecinvert [measure center $ref]]
     }
}
```

スクリプトとして使うべし。

使い方は以下の通り。

`center "all"`

`center "resname POPE"`

### 分子の重心を測る

```
proc center_of_mass {selection} {
	        # some error checking
	        if {[$selection num] <= 0} {
	        error "center_of_mass: needs a selection with atoms"
	        }
	        # set the center of mass to 0
	        set com [veczero]
	        # set the total mass to 0
	        set mass 0
	        # [$selection get {x y z}] returns the coordinates {x y z} 
	        # [$selection get {mass}] returns the masses
	        # so the following says "for each pair of {coordinates} and masses,
	        #  do the computation ..."
	        foreach coord [$selection get {x y z}] m [$selection get mass] {
	           # sum of the masses
	           set mass [expr $mass + $m]
	           # sum up the product of mass and coordinate
	           set com [vecadd $com [vecscale $m $coord]]
	        }
	        # and scale by the inverse of the number of atoms
	        if {$mass == 0} {
	                error "center_of_mass: total mass is zero"
	        }
	        # The "1.0" can't be "1", since otherwise integer division is done
	        return [vecscale [expr 1.0/$mass] $com]
}
```

スクリプトとして使うべし。

使い方は以下の通り。

`center of mass [atomselect top "all"]`

### 作成した構造の保存

`$val writepdb conf.pdb`

conf.pdbは保存するファイル名。確かwritegroとかもできるはず。

### pdbからpsfファイルを作成する

```
package require psfgen
topology ./top_all36_prot.rtf
segment A {pdb input.pdb}
coordpdb input.pdb A
guesscoord
writepdb out.pdb
writepsf out.psf
```

みたいにすればいい。
但しタンパク質の場合はこの限りではない。タンパク質の場合は1本に複数residueがある関係で
重複するresidueを認識できない。
つまりChainAのGLYとChainBのGLYを識別できない。

そのためこの場合は

```
segment A {pdb monomerA.pdb}
coordpdb monomerA.pdb A

segment B {pdb monomerB.pdb}
coordpdb monomerB.pdb B
```

のように各モノマー毎のpdbを用意する必要がある。



## 周期境界関連

### 周期境界を表示する

`pbc box`

### 周期境界を原点周りに表示する

`pbc box -center origin`

### 周期境界のサイズを取得する

`pbc get`

### 周期境界条件を外してunwrapする(主にLAMMPS)

`pbc unwrap`

### 周期境界条件を適用してwrapする(主にLAMMPS)

`pbc wrap -center origin -compound fragment -all`

オプションはお好みで。ちなみにタンパク質は１つの分子内に複数の残基がある都合で

これでは上手くwrapできない。

そのため読み込みに使用するpsfファイルでタンパク質の名前をモノマー単位で変更する。

例えばPRO1,PRO2と言ったように。

その後もう一度VMDを開いて上記のコマンドを実行するとうまくwrapできる。

ちなみにGROMACSを使用しているのであればtrjconvを使う方が早い。



## 便利なコマンド

### GUIを起動せずにVMDを起動する

`VMD -dispdev text out.psf conf.pdb`

### 構造を読み込む

`mol load gro input.gro`

勿論pdbでも可。

### 構造を書きだす

`animate write gro output.gro`

勿論pdbでも可。

これに関して、**atomselect**のオブジェクトに対して使用できる**writepdb**と何が異なるのかと思うかもしれない。

writepdbはオブジェクトに対しての書き出しコマンドであるので、

指定した分子の構造であったりその分子に対して行った変化も記録して反映してくれる。

一方でanimate writeは恐らく読み込んだ元の構造の情報しかない。そのためこのコマンドで書き出しても、

分子を平行移動した結果は構造ファイルに保存されない...はず。
