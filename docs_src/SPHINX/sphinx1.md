# Sphinxの使い方その1

## -Sphinxとは-

### Sphinxって何？Python使えなきゃだめ？

**Sphinx**とはこのウェブページを作成するために使用しているソフトのことである。**Python**を利用して作られたプログラムであるが、ウェブページを作成するためにPythonについて熟知している必要はない。

後述するがSphinxのインストールの際に`pip`を利用する際と、またSphinxのカスタマイズに必要である`conf.py`の記述に使用する程度である。また`conf.py`の記述も自分で考えるというよりは、より便利にするために他の利用者が利用している記述を真似するといった使い方になる。

但しその`conf.py`の記述のために調べることがやや多く、初めて触る人には少し敷居が高いため、ここではSphinxのインストールから実際に自分のページを作るところまでを解説する。またよりSphinxを便利に使うために筆者が設定している`conf.py`や`css`の記述についても記載しておくので参考にしてほしい。

``` caution:: Sphinxは更新の頻度がかなり高い。サポートされているという点ではいいことではあるが、sphinx-quickstartの質問やconf.pyなどの記述がこれまでもわずか2か月ほどで大きく変わってきている。筆者のSphinxのバージョンは3.5.2であるので、バージョンがより新しいものについては適宜調べて利用してほしい。
```

### Sphinxについてもう少し知る(本題)

より具体的には、Sphinxは**reStructuredText**記法で書かれたテキストファイルを**HTML**に変換するソフトである。ウェブページを作成したい場合には、基本的に**HTML**と**CSS**を書く必要がある。

HTMLはテキストボックス、文字、表、画像、リンクなど、WEBサイトに必要な部品を配置するための言語である。またCSSはそのHTMLにデザインを持たせるためものである。要するにかっこいいサイトを作るためにはCSSを書くということである。しかしこれらを習得するには時間がかかる。

それを簡単に行うためのソフトがSphinxであり、そのSphinxにおいてウェブページの元となる***原稿***の書き方がreStructuredText記法(**rst**)である。但し、rstよりももっと簡単にSphinxの原稿を書くための記法があり、それがMarkdown記法(**md**)である。Sphinxはデフォルトではrstにしか対応していないが、**recommonmark**と呼ばれるライブラリをインストールしてconf.pyを設定すれば使用できるようになる(詳細は後述)。

mdあるいはrstで書いた原稿をある種のコンパイルすることによってHTMLとCSSに変換することで、テンプレートはあるものの、かっこいいウェブページが簡単に作れるということである。

もしもっと自由にウェブページを作成したいと思うのであればHTMLやCSSを学ぶことで達成できる。

``` tip:: reStructuredTextやMarkdownのことをマークアップ言語という。
```
### Sphinxを使えなくてもウェブページを作れる

正確にはSphinxを使用できなくても、ウェブページを作るお手伝いができる！と言った方が正しいのかもしれない。先にも述べたようにSphinxを使ったウェブページの作り方は以下のたった2ステップである。

1. mdあるいはrstでページの原稿を書く
2. コンパイルしてmdまたはrstをHTMLに変換

この2番に関しては全員がする必要はなく、Sphinxを使える人にしてもらえばよい。要するに原稿を書いて渡すだけで、その原稿をウェブページに組み込んで貰える。

そのためもしウェブページに興味が少しでもあるのであれば、mdあるいはrstで原稿を書くところから始めてみるといいだろう。



## Sphinxをインストールする

前置きが長くなったが、インストールについて説明していこう。(とはいいつつまた話それるかも...)

Sphinxのインストールのためには**Anaconda**をインストールするのがおすすめ。ちなみにAnacondaとはPythonのディストリビューションの1つであり、ディストリビューションとはコンパイルしてある設定済みのソフトウェアの集合体を指す、らしいが分かりにくい。もう少しかみ砕く。

Pythonでは便利なライブラリが多数用意されているために高度な数値計算や機械学習またMDの解析なども簡単に行える。しかしそれらの多くはPythonに標準で付属しておらず、別途自分でインストールする必要があったりする。その手間を省いてくれるのがAnacondaであり、AnacondaをインストールするだけでPython本体とライブラリがインストールされるため、環境構築がかなり楽になる。Anacondaのインストール(Ubuntuの場合の)方法は以下のサイトが分かりやすく説明しているのでそちらを参考にしてほしい。

[Anacondaのインストール](https://www.pc-koubou.jp/magazine/38846)



実際Anacondaをインストールするとその時点でSphinxがもう利用できるようになっているので、Sphinxのインストールは完了である。

という事実に恥ずかしながら最近気づいた。



というのも筆者はPythonのパッケージを管理するpipコマンドを使用してインストールした。但しpipコマンドがデフォルトではなかったため、わざわざaptを使用してpipを入れた後にpipでSphinxをインストールするという手間を行っていた。まあPythonを本格的に利用するのではなく、Sphinxのためだけに使用するというのであればこれでいい気もする。

この方法でUbuntuにインストールする場合は以下の手順。

```shell
sudo apt install python3-pip
pip3 install sphinx
```

なんにせよ、これでsphinxのインストールができたはずだ。

確認として`sphinx-quickstart`コマンドを実行して

```
Welcome to the Sphinx 3.5.2 quickstart utility.

Please enter values for the following settings (just press Enter to
accept a default value, if one is given in brackets).

Selected root path: .

You have two options for placing the build directory for Sphinx output.
Either, you use a directory "_build" within the root path, or you separate
"source" and "build" directories within the root path.
> Separate source and build directories (y/n) [n]:
```

などと出力されればインストールに成功している。

もしインストールが成功していそうにも関わらず、コマンドが見つからない場合はパスがちゃんと通っているか確認してほしい。sphinxのコマンドは`~/.local/bin`の中にあるはずだ。

## -Sphinxを実際に使ってみる-

上記のように`sphinx-quickstart`コマンドを実行することで、すでにSphinxの利用を開始している。`sphinx-quickstart`により1つのプロジェクトを開始した(1つのホームページ作りを開始した)という状況である。これによりウェブページ作成に必要なディレクトリ群がつくられる。そのため予めウェブページ用のディレクトリを作成してその中で`sphinx-quickstart`とするとよいが、ディレクトリ群がつくられた後でもそれらをまとめてどこかのディレクトリに移せば問題はない。

`sphinx-quickstart`を実行するといくつか質問をされるが、初めての場合どう答えたらよいか迷うと思うので以下に例を載せておく。

``` tip:: sphinxのバージョン2.x.xなどではこの質問の数がとても多く煩わしかったが、3.x.xでは簡単になった。
```

```
Welcome to the Sphinx 3.5.2 quickstart utility.

Please enter values for the following settings (just press Enter to
accept a default value, if one is given in brackets).

Selected root path: .

You have two options for placing the build directory for Sphinx output.
Either, you use a directory "_build" within the root path, or you separate
"source" and "build" directories within the root path.
> Separate source and build directories (y/n) [n]: n

The project name will occur in several places in the built documentation.
> Project name: webpage
> Author name(s): I.K.
> Project release []:

If the documents are to be written in a language other than English,
you can select a language here by its language code. Sphinx will then
translate text that it generates into that language.

For a list of supported codes, see
https://www.sphinx-doc.org/en/master/usage/configuration.html#confval-language.
> Project language [en]: ja
```

何も入力していないところはEnterを押してデフォルトの値を使用している。

質問に答え終わると、ファイルやらディレクトリやらがつくられている。

主要なもののみ簡単に説明しておく。

* Makefile :　Linuxでよく目にするであろうMakefileと同じ。`make html`の実行に必要。
* _build :　　`make`するとこの中にhtmlファイルなどが作られる。
* _static :　　Sphinxで作るページを自分好みの色やデザインにするためにこの中にCSSファイルを作る。
* index.rst :　最も重要。すべてのページの目次であり、ここに他のmdファイルやrstファイルのパスを指定することで自動的にまとめてページの構成を管理してサイドバーや、「次へ」や「前へ」のリンクを作成してくれる。ちなみに拡張子からも分かるようにindex.rst自体もrst記法で書かれたrstファイルである。
* conf.py :　　index.rstの次に重要。このファイルで関数の定義やモジュールのimportをしておくことでSphinxが格段に使いやすくなる。


