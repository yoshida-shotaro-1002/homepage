# Linux commands (nessesary)

先の章で説明したように、linuxは`cd`や`mkdir`などのコマンドを用いて操作します。

このページではその他のコマンドについて解説します。



```eval_rst
=====================================  =====================================
 **command**                             **result**
 pwd                                    カレントディレクトリのパスを表示
 cd                                     ディレクトリ移動
 mkdir                                  新規ディレクトリ作成
 ls                                     ディレクトリの中身表示
 cp                                     ファイル・ディレクトリのコピー
 mv                                     ファイル移動 or ファイル名変更
 rm                                     ファイル・ディレクトリ削除
 echo                                   文字の表示
 cat                                    ファイルの中身を表示
 sl                                     ここまで疲れたらやってみてね
=====================================  =====================================
```

ここでは、

```eval_rst
:doc:`/linux/linux1` で紹介しなかったものについて解説をします。
```



## `cp`

ファイルをコピーするコマンド

### 文法

```shell
cp コピー元ファイル　コピー後名
```



### test1.dat をtest2.datという名前でコピーする。`cp test1.dat test2.dat`

```shell
$ ls
text1.dat

$ cp text1.dat test2.dat # test1.dat をtest2.datという名前でコピーする。

$ ls
text1.dat test2.dat
```



### fileディレクトリの中のtext1.datを一つ上のディレクトリにコピーしてくる。`cp file/text1.dat ./`

```shell
$ ls
file

$ cp file/text1.dat ./

$ ls
file text1.dat
```



### fileディレクトリをfile_sameという名前でコピーする。`cp -r file file_same`

```shell
$ ls
file text1.dat

$ cp -r file file_same

$ ls
file file_same text1.dat
```

このように、ファイルなら`cp`で問題ないがディレクトリをコピーする場合、`cp -r`とする。



## `mv`

- ファイル名を変更するコマンド
- ファイルを移動するコマンド

### ファイル名変更時の文法

```shell
mv 変更前ファイル名　変更後ファイル名
```



### text1.dat をtext_new.datに改名`mv text1.dat text_new.dat`

```shell
$ ls
text1.dat
$ mv text1.dat text_new.dat

$ ls
text_new.dat
```

------



### ファイル移動時の文法

```shell
mv 移動したいファイル　移動先のディレクトリ
```



#### text1.dat をfileの中に入れる`mv text1.dat file`

```shell
$ ls
file text1.dat

$ mv test1.dat file

$ ls
file

$ ls file/
test1.dat
```

移動先のディレクトリはすでに存在している必要がある。まだ作ってなければtext1.datがfileに改名されてしまう。



## `rm`

ファイルを削除するコマンド

### 文法

```shell
rm 削除したいファイル
```



### ファイル test.dat を削除`rm test.dat`

```shell
$ ls
test.dat

$ rm test.dat

$ ls

```

このように、確認メッセージなどがなく、削除されるため注意が必要！！

（ゴミ箱に行ったりもせずいきなり消える）



### ディレクトリ file/ を削除 `rm -r file`

```shell
$ ls
file

$ rm -r file

$ ls

```

`cp` のと同様にディレクトリに対して行うには`-r`のオプションが必要。



## `echo`

文字の表示をするコマンド

### 文法

```shell
echo 表示したい文字
```

###

### I am a geniusと表示する。`echo "I am a genius"`

```shell
$ echo "I am a genius"
I am a genius
```

これだけだと何の面白みもないが、あとで紹介する「リダイレクト」を活用することで利便性が格段に上がる。



### I am a geniusという内容が書かれたtest.datというファイルを作成。`echo "I am a genius" > test.dat`

```shell
$ echo "I am a genius" > test.dat
$ cat test.dat
I am a genius
```

`cat `はファイルの中身を表示するコマンドで、次に説明する。



## `cat`

ファイルの中身を表示する。

### 文法

```shell
cat 中身を表示したいファイル
```



### I am a geniusという内容が書かれたtest.datというファイルの中身を見る`cat test.dat`

```shell
$ cat test.dat
I am a genius
```



### ~/.bashrcの中身を表示してみる。

```shell
$ cat ~/.bashrc
........................
........................
```



## `sl`

ここまでやって疲れたら`sl`コマンドを試してほしい。疲れがなくなるはずだ。

