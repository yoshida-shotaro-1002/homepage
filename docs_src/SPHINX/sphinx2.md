# Sphinxの使い方その2

## -ウェブページの作成-

以下のコマンドを、**Makefile**と同じ階層のディレクトリで実行してみてほしい。

```shell
make html
```

その後windowsのエクスプローラーで_build/htmlの中にある**index.html**をダブルクリックしてほしい。

すると以下のようなサイトがブラウザで見られるはずだ。


```eval_rst
.. figure:: image/first_indexrst.png
```

ここまでの手順で何もページ作りをしていないが、デフォルトのページが作れたことになる。

何も作っていないのに何を基にHTMLを作っているのかと疑問に思う人もいるかもしれないが、前ページの[「その1」](/SPHINX/sphinx1.md)で述べたようにindex.rstが目次の役割を果たしつつ、それ自体もrstの原稿であったことを思い出してほしい。したがってindex.htmlはindex.rstをコンパイルして作られたページということになる。実際index.rstの中身は以下のようになっており、

```
.. webpage documentation master file, created by
   sphinx-quickstart on Thu Apr 29 02:53:00 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to webpage's documentation!
===================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
```

Welcome to webpage's documentation!という見出しの言葉が対応していることがわかるだろう。

``` tip:: rstでは「=」を見出しにしたい文字列の下に連続して書くことで、その文字列を見出しにできる。
```

```eval_rst 
.. important::
    また.. hogehoge:: はSphinx特有の書き方である。中でも上記のtoctreeはSphinxの大黒柱ともいうべき機能であり、これによりページ間の繋がりを管理して目次等を作ってくれている。
```

余談だが、索引をクリックすると索引と呼ばれるページが出てくるもののこのページの使い方を筆者は理解できていない。というか必要ないように思える。またモジュール索引に関してはクリックするとリンクが切れており、別途設定が必要なのであろう。さらに検索ページも、そもそも開いたページに検索フォームが作られているため必要ない。

したがって筆者はこれらのページを全て削除している。

すなわち`Indices and tables`以下の行を全てindex.rstから消している。以降の説明でも同様であるので注意してほしい。

## -実際に自分で書いてみる-

Hello Worldという見出しが表示されるだけのページを作成してみよう。

まずはrstファイルを作る。作る場所はindex.rstと同じ階層に作成する。

すなわち筆者の環境では以下のようなディレクトリ構造になっている。

```
.
├── Makefile
├── _build
├── _static
├── _templates
├── conf.py
├── index.rst
├── hello.rst
└── make.bat
```

rstはただのテキストファイルであるので、rstファイルを書くためのエディタはviでもメモ帳でもVScodeでも何でもよい。以下のように記述する。

```
.. hello.rst

Hello World
============
```

ちなみに`.. hello.rst`はコメントであり、Sphinxではこのようにしてコメントアウトできる。

``` tip:: rstのコメントアウトの仕方ではなく、あくまでSphinxのコメントアウトである。ただrstの記法はそもそもSphinxでしか用いられていなかった気もする。例えばmd記法はSphinxで使える以外にも様々なところで使える。例えばqiitaではmd記法を使用してページを作成できる。
```

hello.rstを作成したら、index.rstにhello.rstのパスを記述する必要がある。作成したいかなるrstやmdもindex.rstにそのパスを記述する必要がある。そうでないとページを管理できないため。

したがってindex.rstに以下のように記述。

```
.. webpage documentation master file, created by
   sphinx-quickstart on Thu Apr 29 02:53:00 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to webpage's documentation!
===================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   hello.rst
```

ちなみにtoctreeにファイルのパスを記述する際はインデントを必ず意識してほしい。今は、toctreeの「t」と同じ位置にhello.rstの「h」が来ている。もしhello.rstの前にスペースあるいはタブを入れずに「hello.rst」と書くと、うまく読み込めないためにエラーが表示される。

ここまでできたら`make html`を実行する。そしてindex.htmlをブラウザで開くと以下のようになっているはずだ。

```eval_rst
.. figure:: image/helloworld1.png
```

```eval_rst
.. figure:: image/helloworld2.png
```

ただしrstファイルが増えていくとindex.rstと同じ階層のディレクトリがめちゃくちゃになる。その場合は、ディレクトリを作ってディレクトリも含めたパスをtoctreeに記述してやればよい。次のようなディレクトリ構図を想定する。


```
.
├── Makefile
├── _build
├── _static
├── _templates
├── conf.py
├── index.rst
├── Hellow
│　　　└──hellow.rst
└── make.bat
```
この場合toctreeには次のように書く。

```
.. webpage documentation master file, created by
   sphinx-quickstart on Thu Apr 29 02:53:00 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to webpage's documentation!
===================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   Hello/hello.rst
```

これができれば一応どんな記事も書けるということになる。
次のページではSphinxを使用したページデザインの変更について説明していく。

``` warning:: rstファイルやmdファイルには必ず見出しを含む必要があり、見出しのないただのテキストだけのファイルは許されずエラーとなる。
```
