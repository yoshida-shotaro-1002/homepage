# Sphinxの使い方その4

## -Markdownを使おう-

さて、何度も言っているようにMarkdownをSphinxで使えるようにしたい。但し正直Sphinxを使用していてreStructureTextもそこまで難しくない気もしており筆者的にはどちらでもいい気もしている。但しマークアップ言語の主流はやはりMarkdownであり、どうせならmdで書いた方がいい気もする。

ということでMarkdownを使えるようにしていく。

SphinxでMarkdownを使用するためにはrecommonmarkというライブラリをインストールする必要がある。以下のコマンドでインストールが可能。

```
pip3 install recommonmark
```

インストール終了後、次のようにconf.pyの`extensions = []`となっているところを`extensions = ['recommonmark']`にする。

```python
# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
        'recommonmark'
]
```

さらに次のように`source_suffix`を記述する。

```python
# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
#html_theme = 'alabaster'
html_theme = 'sphinx_rtd_theme'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']
html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]

source_suffix = {
        '.rst': 'restructuredtext',
        '.txt': 'markdown',
        '.md' : 'markdown',
}
```

これはファイルの拡張子がrstの時、txtの時、mdの時、それぞれどういうファイルとして読み込みを行うかを記述している。この記述によりmd記法だけではなく元々のrst記法も使用できるようにしている。

さらにその下に次のような記述をする。

```
from recommonmark.transform import AutoStructify

github_doc_root = 'https://github.com/rtfd/recommonmark/tree/master/doc/'
def setup(app):
        app.add_config_value('recommonmark_config', {
                'url_resolver': lambda url: github_doc_root + url,
                'auto_toc_tree_section': 'Contents',
                }, True)
        app.add_transform(AutoStructify)
```

正直筆者はこの内容を把握していないが、recommonmarkを使用するためのメインの記述という理解で十分だと思う。一応参考にしたサイトのリンクを記載しておく

[Sphinxでmarkdownを扱う(qiita)](https://qiita.com/leo-mon/items/46c43f0f97f730e64754)

ただし上のリンクで`from recommonmark.parser import CommonMarkParser`などとしているが、Sphinxのバージョンが3.x.xになって以降はparserに関する記述は必要なくなったはず。そのため上の説明ではこの部分の記述はしていない。逆に古いSphinxを使用している場合は必要に応じて調べてほしい。

また`setup()`の意味が知りたければrecommonmarkの公式ドキュメントを見るとよい。詳しく書いてあるはずだ。

[recommonmark公式ドキュメント](https://recommonmark.readthedocs.io/en/latest/index.html#autostructify)

ただしこの辺りの記述はSphinx特有のものではなくゴリゴリにPythonなのでPythonを軽く学んだ後に読むのがよいかもしれない。

話がそれたがこれでmd記法を使用するための準備は完了である。

mdでファイルを作成し、index.rstにパスを記述した上で`make html`を実行してちゃんとコンパイルができるかどうか確認してみてほしい。

参考までにmdファイルとindex.rstの中身を記載しておく。

```
[](
hello2.md
)

# Hello World
```

hello2.mdの記述はコメントである。またmdとrstで拡張子より前の部分が同じ名前のファイルが複数あると警告が出るため、hello2.mdとした。

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
   Hello/hello2.md
```

これでMarkdownがSphinxで使用できるようになった。

## -その他のカスタマイズ-

### CSSで見栄えをよくする

もう十分設定したしお腹いっぱいな気もする。だがここまできたら次で説明する設定もしてほしい。ほかに何を設定するかというと、例えば自分で作ったサイトをもう１度見てもらうと分かる。[「その3」](/SPHINX/sphinx3.md)のrtdのテーマの写真をもう1度載せておく。

```eval_rst
.. figure:: image/theme_rtd.png
```

上の画像を見ると、細かいことかもしれないがウェブページの端がパソコンの端に一致していないことが分かる。よく分からないと思った人は、ノートPCではなくデスクトップで見てほしい。なぜそんなにページが縮こまっているんだ...と思うはず。これはrtdのテーマでサイトの端が固定されていることが原因。

端の固定を解決するには**CSS**を記述し、そのCSSを反映するような形でSphinxでコンパイルしてもらえばいい。

とは言ってもCSSを1から書くのではなく直したいところのみを記述すればいい。CSSを上書きするというイメージの方が正確なのかもしれない。

そのためにはまず、`_static`ディレクトリにcustom.cssというファイルを作る。さらにファイルに次のような記述をする。

```css
.wy-nav-content {
        max-width: none;
}
```

この記述はCSS特有の記述方法である。

さらにMarkdownを使用するために定義した`setup()`にCSSに関する内容も書いてやる。つまりconf.pyは以下のようになる。

```python
from recommonmark.transform import AutoStructify

github_doc_root = 'https://github.com/rtfd/recommonmark/tree/master/doc/'
def setup(app):
        app.add_config_value('recommonmark_config', {
                'url_resolver': lambda url: github_doc_root + url,
                'auto_toc_tree_section': 'Contents',
                }, True)
        app.add_transform(AutoStructify)
        app.add_css_file('custom.css')
```

こうすることでページを作る際にcustom.cssも参照してページが作られるようになる。

`make html`を実行しよう。下の画像のように横幅がちゃんと広がっているはずだ。

```eval_rst
.. figure:: image/fix_rtd.png
```

もう１つCSSのカスタマイズを説明しておく。次のような記述をcustom.css足してみる。

```
.wy-side-nav-search, .wy-nav-top {
        background: #32CD32
}
```

その後`make html`を実行する。

```eval_rst
.. figure:: image/green_rtd.png
```

見ての通りページのメインカラーが変わる。地味な変化かもしれないが、応用すればサイドバーや文字など様々なところのデザインを自分好みに変えることができるはずである。また検索ボックスの上に現在は家のマークとwebpageという文字が書いてあるが、この部分を好きなロゴに変えることができるらしい。T先輩が変えていたので間違いない。そのためCSSやSphinxに関して色々調べて触ってみるとよい。

``` hint:: #32CD32はCSSにおいて緑色に割り当てられている番号である。色の種類は番号の数だけあり、細かく変更できる。「CSS　色　番号」などでググれば色々出てくると思うので調べて自分好みの色にするとよい。
```

### Markdownの表組み

このページトップでrecommonmarkによりMarkdownを使用できるようにした。しかしながらrecommonmarkはMarkdownの表組み記法に対応していないらしい。そのためMarkdownで表を書きたい場合は別途ライブラリが必要になる。インストールは次のコマンド。

```
pip3 install sphinx-markdown-tables
```

インストールが終わったら、conf.pyに`sphinx-markdown-tables`を利用することを明示する。やり方は`extensions`に`sphinx_markdown_tables`の記述を足す。つまり以下のようにすればいい。

```
extensions = [
        'recommonmark',
        'sphinx_markdown_tables'
]
```

これで設定は完了している。

### 著者・プロジェクト名など

はじめの質問で著者やプロジェクト名を考えるのに悩んだ人もいるかもしれない。しかしこれらは後からも変更ができるので安心してほしい。変え方は`project`,`copyright`,`author`の変数にそれぞれ文字列を代入するだけである。

```
project = 'webpage'
copyright = '2021, I.K.'
author = 'I.K.'
```



***



ここまで色々設定してきたため、conf.pyの中身がよく分からなくなってしまったという人もいるかもしれない。そのため筆者が現在使用しているconf.pyの中身を記載しておく。

```
# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))
import sphinx_rtd_theme


# -- Project information -----------------------------------------------------

project = 'MDnotes'
copyright = '2021, S.Harada,I.Kawabata,S.Ikeda,K.Shibata'
author = 'ikawabata'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
        'recommonmark',
        'sphinx_markdown_tables'
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# The language for content autogenerated by Sphinx. Refer to documentation
# for a list of supported languages.
#
# This is also used if you do content translation via gettext catalogs.
# Usually you set "language" from the command line for these cases.
language = 'ja'

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
#html_theme = 'alabaster'
html_theme = 'sphinx_rtd_theme'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']
html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]

source_suffix = {
        '.rst': 'restructuredtext',
        '.txt': 'markdown',
        '.md' : 'markdown',
}

from recommonmark.transform import AutoStructify

github_doc_root = 'https://github.com/rtfd/recommonmark/tree/master/doc/'
def setup(app):
        app.add_config_value('recommonmark_config', {
                'url_resolver': lambda url: github_doc_root + url,
                'auto_toc_tree_section': 'Contents',
                }, True)
        app.add_transform(AutoStructify)
        app.add_css_file('custom.css')
```

設定についてはここまで。次のページからは高度なmdとrstの書き方について説明する。
