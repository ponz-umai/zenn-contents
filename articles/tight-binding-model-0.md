---
title: "Tight-binding modelの第二量子化表示を初歩からちゃんと理解する"
emoji: "👶"
type: "tech" # tech: 技術記事 / idea: アイデア
topics: ["quantum","quantumcomputing","quantumcomputer","物理","物理学"]
published: true
---
# Tight-bindingについて
物性物理では以下のような「飛び移り積分」（や"Hopping Integral", "Transfer Integral"等）と呼ばれる謎の定数$t$と、サイト$i,j$に電子を作ったり消したりする謎の演算子「生成・消滅演算子」とその突然のエルミート共役とあと謎の和の取り方$\sum_{\left<ij\right>}$等によって構成される謎の式いわゆる第二量子化表示のTight-binding modelがよく使われます：

$$
\mathcal{H} = \sum_{\left <i,j \right>,\sigma }-t_{ij}\left( \hat{a}^\dagger_{i\sigma }\hat{a}_{j\sigma }   + h.c.\right).
$$

このモデルは固体内の電子が格子点（原子核や分子イオン）に強く束縛され、ほぼ局在していると仮定した描像に対応しているらしいです。

Tight-binding modelは非常にシンプルな構造ですが、そこに電子間のオンサイトクーロン相互作用を加えたHubbard model

$$
\mathcal{H} = \sum_{\left <i,j \right>,\sigma }-t_{ij}\left( \hat{a}^\dagger_{i\sigma }\hat{a}_{j\sigma } + h.c. \right)  + U\sum_i\hat{n}_{i\uparrow}\hat{n}_{i\downarrow}
$$

に始まり様々なモデルの出発点となっており、それらからMot絶縁体やBCS超伝導、Anderson局在、はたまた今ホットなトポロジカル物性、量子ビット候補のマヨラナ粒子、、等々様々な物理現象が描き出されることになり、物性物理の勉強・研究をしていれば必ず出会うことになります。または物性とは異なる分野の方でも目にする機会はあるのではないかと思います。Hubbard modelの面白さについては例えばこちらの[田崎先生のpdf](https://www.gakushuin.ac.jp/~881791/pdf/KBHubbard.pdf)等を読むとわかりやすく幅広く書いてあります。


# っていうことなのですが
ひとたび上記のような生成・消滅演算子で記述されたモデル（ハミルトニアン）を受け入れてしまえば、その先は演算子の交換関係や、平均場近似等々を使ってある程度機械的にモデルの解析を進めていくことができます。

しかし「サイト$i,j$の電子って何？」とか、「飛び移り積分って？」とか「電子を作ったり消したりする演算子っていったい何？？？」とか考えてみるとなんだかよくわからない、（学部時代を適当に過ごしてきた）大学院生や、他分野の研究者・技術者の方もいるのではないでしょうか。
というか何を隠そう、私が学部時代に量子力学を適当に勉強した状態で修士に進み、その辺にもやもやを抱えながら研究してきた残念院生だったのですが（そしてそれを解消できないまま修了しちゃったのですが）、このままだと死んでも死にきれないので過去の自分がブラックボックスにしてきた部分を理解するべく勉強した内容をまとめようと思い立ちました。


量子力学についてのなんとなくの理解を前提に、上記Tight-binding modelの第二量子化表示について納得できるまでの内容を順次作成していこうと思います。とはいえ全部書くとすごい量になってしまい終わらなさそうなので、適宜既存の教科書に説明を譲りつつ、教科書では省略されていたりあまり触れられていない部分を中心に書いていこうと思っています。

世の中には「第二量子化表示の説明で終わる量子力学の教科書」と、「第二量子化には触れずにTight-binding近似について説明してある物性の教科書」と、「第二量子化は当然知っているものとして突然上記のモデルが登場する物性の教科書」の3種類しかない気がするので、そのあたりのギャップを埋められたらと思っています。（と修士の頃は思っていたのですが、最近改めて探すと色々と親切な教科書が見つかったりもするのですが）

→2023年2月10日　ようやく完結しました。長かった。。。手っ取り早く結論が知りたい方は[Tight-binding modelの第二量子化表示](https://zenn.dev/ponzumai/articles/tight-binding-model-2nd-q-tbmodel)の章をご覧ください。一つ一つ進めていきたい方は以下の目次の順番で見ていっていただければと思います。

# 目次（それぞれ一つの記事のリンクになってます。）


## [水素原子中の電子](https://zenn.dev/ponzumai/articles/tight-binding-model-hydrogen-atom)
## [（補足）水素原子中の電子のシュレディンガー方程式の具体的な計算](https://zenn.dev/ponzumai/articles/appendix-hidrogen-atom)
## [多電子状態の波動関数](https://zenn.dev/ponzumai/articles/tight-binding-model-many-electron)
## [電子のスピンを考慮した多電子系の波動関数](https://zenn.dev/ponzumai/articles/tight-binding-model-spin)
## [多電子原子中の電子](https://zenn.dev/ponzumai/articles/tight-binding-model-many-electron-atom)
## [格子ベクトルと逆格子ベクトル](https://zenn.dev/ponzumai/articles/tight-binding-model-lattice-vec)
## [逆格子ベクトル・逆格子点の物理的な意味](https://zenn.dev/ponzumai/articles/tight-binding-model-reciprocal-lattice)
## [固体（結晶）中の電子状態](https://zenn.dev/ponzumai/articles/tight-binding-model-electrons-in-solids)
## [Bloch関数の局在関数ーWannier関数ーを用いた展開](https://zenn.dev/ponzumai/articles/tight-binding-model-wannier-func)
## [第一量子化のTight-bindingモデル（前編）](https://zenn.dev/ponzumai/articles/tight-binding-model-1st-q-1)
## [第一量子化のTight-bindingモデル（後編）](https://zenn.dev/ponzumai/articles/tight-binding-model-1st-q-2)
## [飛び移り積分（Hopping Integral）の物理的意味・Wannier関数の従う方程式](https://zenn.dev/ponzumai/articles/tight-binding-model-hopping-int)
## [重なり積分（Overlap Integral）の物理的意味（の妄想）](https://zenn.dev/ponzumai/articles/tight-binding-model-overlap-int)
## [第二量子化表示/生成・消滅演算子の方法](https://zenn.dev/ponzumai/articles/tight-binding-model-2nd-q)
## [Tight-binding modelの第二量子化表示](https://zenn.dev/ponzumai/articles/tight-binding-model-2nd-q-tbmodel)


# さいごに
質問やわかりにくい点などあればコメントを頂けますと幸いです。（もちろん感想だけでもいただけるととても喜びます。）
また、上述の通り滑り込み修士が独学をベースに書いているノートなので、勘違いや間違いは多数含まれている可能性がありますし、おかしな書き方をしている部分もあるかと思いますのでご注意ください（指摘していただければ嬉しいです）。