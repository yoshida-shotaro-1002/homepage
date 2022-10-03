# 高分子鎖の物性値調整に対する準備

## 物性値の調整に必要なプログラムを用意する

NAS上においてあるので/mnt/Impact02/temp_yoshidaのディレクトリに入ってAdjustEntanglement(からみあい調整用のディレクトリ)、AdjustRg(慣性半径調整用のディレクトリ)、bin(物性値調整の実行ファイルなどがまとめられているディレクトリ)をコピーして持ってくる。

## 上記のプログラムを動かすのに必要なアプリを用意する。

先ほどのbinの実行ファイルを動かすためにnetcdfとliblbfgsが必要である。手元でプログラムを動かす分には他のページでこれらをインストールしてPATHを通してあると思うので問題ないが、スパコンでプログラムを動かす過程が存在するためスパコンにもnetcdfとliblbfgsを入れてPATHを通しておいて欲しい。このやり方は、基本的に手元のパソコンにインストールする方法と同じであるので割愛する。補足として解凍前の〜.tar.gzファイルをスパコンにコピーする際に以下のようにrsyncコマンドでアップロードする。

```shell
# LIBLBFGSを手元のパソコンからスパコン(分子研)にアップロードする場合
rsync -avh ~/apl/liblbfgs-1.10.tar.gz RCCS:/home/users/dfw/2022/apl/
# dfwは筆者のスパコン上でのユーザー名である。RCCSは筆者の手元のパソコンにおけるconfigでの分子研にあるスパコンの登録名
```

また、筆者のスパコン上の.bashrcファイルは以下のようにPATHを通してある。(該当部分のみ掲載)

```shell
# User specific aliases and functions
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/users/dfw/2022/apl/liblbfgs-1.10/install/lib/:/home/users/dfw/2022/apl/netcdf-3.6.3/install/lib/:/local/apl/lx/gromacs2018/lib64/
```
