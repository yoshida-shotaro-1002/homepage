# Sphinxの使い方その3

## -Sphinxのデザインを変更してみる-

Sphinxにはページデザインのテーマが存在し、前のページで作成したウェブページも1つのテーマを使って作成している。しかしSphinxには他にもテーマが組み込まれている。テーマには次のようなものがある。

- alabaster
- default
- sphinxdoc
- scrolls
- agogo
- nature
- pyramid
- haiku
- traditional
- bizstyle

テーマの変更はconf.pyを編集することで変更できる。conf.pyでテーマに関する記述は以下のようになっている。

```
# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'alabaster'
```

現在はalabasterに設定されているためこの変数を上記のテーマの名前に変更すればよい。例えばdefaultに変更してみるとページは次のようになる。

```eval_rst
.. figure:: image/theme_default.png
```

これは**GROMACS**の公式ドキュメントのデザインと同じものである。またsphinxdocの場合は次のようになる。

```eval_rst
.. figure:: image/theme_sphinxdoc.png
```

他のテーマについても色々試してみてほしい。

## -Read the Docs theme for Sphinx-

Sphinxには様々なデザインがあるし色々紹介しておいて言うことではないかもしれないが、どれも少し古いデザインな気もする。そんな人には**Read the Docs theme for Sphinx**の利用をお勧めする。これはSphinxではなく第三者がSphinxで使えるように製作したSphinxのテーマの1つである。正直かっこいいし、T先輩がこのデザインを使用していたというのもある。但し先ほどまでのようにconf.pyの変数を変えるだけでは利用できないので、準備の仕方を説明していく。

まずはRead the Docs theme for Sphinxをインストールする。コマンドはこちら。

```
pip3 install sphinx_rtd_theme
```

ちなみにAnacondaをインストールすることでSphinxをインストールしたから、[「その1」](/SPHINX/sphinx1.md)でpipなんてコマンド使わなかったよという人も安心してほしい。Anacondaの中にちゃんとpipコマンドも入っているので、Anacondaにパスが通っていればちゃんと使用できるはずだ。

その後次の記述をconf.pyに足してやればよい。

```
import sphinx_rtd_theme
html_theme = "sphinx_rtd_theme"
html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]
```

conf.pyのどこに書いても大丈夫だと思うが、筆者は次の位置に記述している。

```python
# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))
import sphinx_rtd_theme
```

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
```

以上のことができたら`make html`を実行してみよう。

次のようなページに変わっているはずだ。

```eval_rst
.. figure:: image/theme_rtd.png
```

ここまでできればSphinxを既に十分活用できていると言えるだろう。ただし、まだまだカスタマイズできる部分はある。次のページではそれを説明していく。

また初めにMarkdownで書けるから楽と言ったがまだこの説明では1度も書いていない。それはconf.pyの設定が少し複雑であったためである。しかし先ほどのrtdテーマのカスタマイズでconf.pyの中身を見ている読者ならそれほど抵抗なくできるはずだ。Markdownを利用するための方法も次のページで解説していく。
