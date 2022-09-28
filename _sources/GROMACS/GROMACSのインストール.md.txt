# GROMACSのインストール

ここではGROMACSのインストール方法をまとめている。 
``` important:: GROMACSのバージョンはgromacs-XXXとして表記しているので、基本的には最新のものを選んでインストールする。
```

## インストール手順

まず、GROMACSをダウンロードする必要がある。 [GROMACS documentation](https://manual.gromacs.org/documentation/)にアクセスし、適切なバージョンのDownloadから「gromacs-XXX.tar.gz」をダウンロードする。 または`wget`を用いてターミナルで直接ダウンロードしてもよい。

```shell
# GROMACSのインストール先へ移動
# インストール先のディレクトリにgromacs-XXX.tar.gzを用意しておく
cd ~/App
# ターミナルで直接ダウンロードする場合はここで
# wget  https://ftp.gromacs.org/gromacs/gromacs-XXX.tar.gz
tar -xzvf gromacs-XXX.tar.gz
cd gromacs-XXX
```

gromacs-XXXのディレクトリへ移動したら`ls`で中身を確認しておく。

次に、ビルド用のディレクトリとGROMACSをcmakeするディレクトリを作成する。

```shell
mkdir build
mkdir gmx_install
```

作成したbuildディレクトリへ移動し、`cmake`を行う。 この時、GROMACSの様々なオプションを指定してインストールすることが多い。 ここではその一例を示す。

| GROMACSインストールオプション |                            説明                            |
| :---------------------------: | :--------------------------------------------------------: |
|    -DGMX_BUILD_OWN_FFTW=ON    | GROMACSがFFTWを自動的にダウンロード、ビルドすることを許可  |
|   -DREGRESSIONTEST_DOWNLOAD   |   GROMACSのテストをダウンロードする（make checkで実行）    |
|  -DCMAKE_INSTALL_PREFIX=xxx   |                GROMACSをcmakeするパスの指定                |
|         -DGMX_MPI=ON          |                並列計算を利用してビルドする                |
|         -DGMX_GPU=ON          |              GPU版のGROMACSをインストールする              |
|    -DCMAKE_C_COMPILER=xxx     |  使用するCコンパイラのパスの指定（環境変数CCでも設定可）   |
|   -DCMAKE_CXX_COMPILER=xxx    | 使用するC++コンパイラのパスの指定（環境変数CXXでも設定可） |
|        -DGMX_DOUBLE=ON        |                 GROMACSを倍精度で構築する                  |
|        -DGMX_SIMD=xxx         |         GROMACSを実行するノードのSIMDレベルの指定          |

ローカルのPCにインストールする場合は系の構築や解析に使用する程度であり、大きな計算には使用しないので、次のオプションを指定すればよいと思う。

```shell
# ローカルでのcmakeの例
cd builds
cmake ../ \
-DGMX_BUILD_OWN_FFTW=ON \
-DREGRESSIONTEST_DOWNLOAD=ON \
-DCMAKE_INSTALL_PREFIX=${HOME}/App/gromacs-XXX/gmx_install \
-DGMX_MPI=OFF \
-DGMX_DOUBLE=OFF \
-DGMX_GPU=OFF
```

一方、スパコンなどにインストールする場合は並列計算を使用するため`-DGMX_MPI=ON`とする。GPUが利用できる場合は`-DGMX_GPU=ON`も指定しておく。

```shell
# スパコンでのcmakeの例
cd build
cmake ../ \
-DGMX_BUILD_OWN_FFTW=ON \
-DREGRESSIONTEST_DOWNLOAD=ON \
-DCMAKE_INSTALL_PREFIX=${HOME}/App/gromacs-XXX/gmx_install \
-DGMX_MPI=ON \
-DCMAKE_C_COMPILER=mpicc \
-DCMAKE_CXX_COMPILER=mpicxx \
-DGMX_GPU=ON \
-DGMX_SIMD=AVX_256
```

cmakeが終了したら続けて`make`を行う。 マルチプロセッサ環境であれば、`-j`を付けて並列コンパイルをすることで高速化ができる。 GROMACSのテストをダウンロードしている場合は同時に`make check`を行う。

```shell
make -j 4
# GROMACSのテストをダウンロードした場合はチェックを実行
make check
```

最後に`make install`をする。

```shell
make install
```

終了するとgmx_installにGROMACSの実行ファイルが生成しているはずである。 .bashrcにパスを通してインストールは完了である。

```shell
.bashrc
# ホームディレクトリの.bashrcを編集
export PATH=${HOME}/App/gromacs-XXX/gmx_install/bin:${PATH}
```

正しくインストールされていれば、次のコマンドでGROMACSのバージョンが表示されるだろう。

```shell
# .bashrcの編集後にパスを更新
source ~/.bashrc
# GROMACSが正しく動作するか確認（バージョン確認）
gmx --version
```
