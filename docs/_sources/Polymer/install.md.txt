# 高分子の研究に必要なライブラリのインストール

ここでは高分子グループにおける研究を進める上で必要となってくるライブラリのインストールについてまとめている。 

## Liblbfgs

まず、Liblbfgsをダウンロードする必要がある。 [Liblbfgs](https://src.fedoraproject.org/repo/pkgs/liblbfgs/liblbfgs-1.10.tar.gz/2a46da6c4014d6b1e8a8790a93edffbb/)にアクセスし、「liblbfgs-1.10.tar.gz」をダウンロードする。 

```shell
# Liblbfgsのインストール先へ移動
# インストール先のディレクトリにliblbfgs-1.10.tar.gzを用意しておく
cd ~/apl
tar -xzvf liblbfgs-1.10.tar.gz
cd liblbfgs-1.10
```

liblbfgs-1.10のディレクトリへ移動したら`ls`で中身を確認しておく。

次に、configure(アプリがインストールされる環境を自動的に調べ、その環境に合わせたMakeFileを自動的に作ること)するディレクトリを作成する。

```shell
mkdir install
```

作成したinstallディレクトリへ移動し、`configure`を行う。 このとき指定がなければ/usr/local内にインストールされる。その際にスーパーユーザーの権限が必要になるため、この後の一連のコマンド入力の頭にsudoをつける必要がある。ローカル内にインストールしたい場合は以下のように--prefixでインストール先を指定する。

```shell
# ローカルでのconfigureの例
./configure --prefix=/home/(ユーザー名)/apl/liblbfgs-1.10/install/
```

configureが終了したら続けて`make`を行う。 マルチプロセッサ環境であれば、`-j`を付けて並列コンパイルをすることで高速化ができる。 

```shell
make -j 4
```

最後に`make install`をする。

```shell
make install
```

終了するとinstallにlibというディレクトリが生成しているはずである。このlibの中にLiblbfgsの実行ファイルが生成している。 .bashrcにパスを通してインストールは完了である。

```shell
vi .bashrc
# ホームディレクトリの.bashrcを編集
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/{ユーザー名}/apl/liblbfgs-1.10/install/lib/
```

## Netcdf

このNetcdfのインストール方法も基本的にLiblbfgsと同様である。まず、Netcdfをダウンロードする必要がある。 [Netcdf](https://src.fedoraproject.org/repo/pkgs/netcdf/netcdf-3.6.3.tar.gz/334e9bdc010b6cd03fd6531a45fe50ad/)にアクセスし、「netcdf-3.6.3.tar.gz」をダウンロードする。 

```shell
# Netcdfのインストール先へ移動
# インストール先のディレクトリにnetcdf-3.6.3.tar.gzを用意しておく
cd ~/apl
tar -xzvf netcdf-3.6.3.tar.gz
cd netcdf-3.6.3
```

netcdf-3.6.3のディレクトリへ移動したら`ls`で中身を確認しておく。

次に、configure(アプリがインストールされる環境を自動的に調べ、その環境に合わせたMakeFileを自動的に作ること)するディレクトリを作成する。

```shell
mkdir install
```

作成したinstallディレクトリへ移動し、`configure`を行う。 このとき指定がなければ/usr/local内にインストールされる。その際にスーパーユーザーの権限が必要になるため、この後の一連のコマンド入力の頭にsudoをつける必要がある。ローカル内にインストールしたい場合は以下のように--prefixでインストール先を指定する。

```shell
# ローカルでのconfigureの例
./configure --prefix=/home/(ユーザー名)/apl/netcdf-3.6.3/install/
```

configureが終了したら続けて`make`を行う。 マルチプロセッサ環境であれば、`-j`を付けて並列コンパイルをすることで高速化ができる。 

```shell
make -j 4
```

最後に`make install`をする。

```shell
make install
```

終了するとinstallにlibというディレクトリが生成しているはずである。このlibの中にNetcdfの実行ファイルが生成している。 .bashrcにパスを通してインストールは完了である。

```shell
vi .bashrc
# ホームディレクトリの.bashrcを編集
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/{ユーザー名}/apl/netcdf-3.6.3/install/lib/
```
## Intel oneapi(basetoolkit, hpctoolkit)

まず、Intelのアカウントを登録する必要がある。[Intel webpage](https://www.intel.com/content/www/us/en/homepage.html/)にアクセスし、右上の人のマークをクリックし、指示に従って登録を進めていく。メールによる認証に時間がかかるため注意する。

Intelのアカウントを作成し、ログインしたらサイト内検索でoneapiと検索する。Base ToolkitとHPC Toolkitが見つかると思うのでLinuxでダウンロードする。

まず、インストールするためのディレクトリを/apl/intel/oneapi/のディレクトリを作成する。l_BaseKit_p_2021.1.0.2659_offline.sh, l_HPCKit_p_2021.1.0.2684_offline.shがダウンロードされているので、それぞれ以下のコマンドを打ち込むことでウィンドウの指示が表示されると思うのが、その指示に従ってインストールを進めていく。

```shell
#Base ToolKitの場合
bash l_BaseKit_p_2021.1.0.2659_offline.sh
#HPC ToolKitの場合
bash l_HPCKit_p_2021.1.0.2684_offline.sh
```
インストールはカスタムしなくても良い。ディレクトリの指示があれば/apl/intel/oneapi/のディレクトリを指定して、両方ともそこにインストールする。
