---
title: "Tight-binding modelの第二量子化表示を初歩からちゃんと理解する"
emoji: "🐈"
type: "tech" # tech: 技術記事 / idea: アイデア
topics: ["quantum","quantumcomputing","quantumcomputer","物理","物理学"]
published: true
---
# Tight-bindingについて
物性物理では以下のようなtTight-binding modelがよく使われます。

$$
\mathcal{H} = \sum_{\left <i,j \right>,\sigma }\left( -t_{ij}a^\dagger_{i\sigma }a_{j\sigma }  \right) 
$$

ここで$\left<i,j \right>$は最近接の格子点について和を取ることを表します。$t_{ij}$は飛び移り積分、$a^\dagger_{i\sigma },a_{j\sigma }$はそれぞれ格子点$i(j)$、スピン$\sigma$の電子の生成（消滅）演算子です。

このモデルは固体内の電子が格子点（原子核や分子イオン）に強く束縛され、ほぼ局在していると仮定した描像に対応しています。

Tight-binding modelは非常にシンプルな構造ですが、そこに電子間のオンサイトクーロン相互作用を加えたHubbard model

$$
\mathcal{H} = \sum_{\left <i,j \right>,\sigma }\left( -t_{ij}a^\dagger_{i\sigma }a_{j\sigma }  \right)  + U\sum_in_{i\uparrow}n_{i\downarrow}
$$

に始まり様々なモデルの出発点となっており、それらからMot絶縁体やBCS超伝導、Anderson局在、はたまた今ホットなトポロジカル物性、量子ビット候補のマヨラナ粒子、、等々様々な物理現象が描き出されることになり、物性物理の勉強・研究をしていれば必ず出会うことになります。または物性とは異なる分野の方でも目にする機会はあるのではないかと思います。Hubbard modelの面白さについては例えばこちらの[田崎先生のpdf](https://www.gakushuin.ac.jp/~881791/pdf/KBHubbard.pdf)等を読むとわかりやすく幅広く書いてあります。


# っていうことなのですが
ひとたび上記のような生成・消滅演算子で記述されたモデル（ハミルトニアン）を受け入れてしまえば、その先は演算子の交換関係や、平均場近似等々を使ってある程度機械的にモデルの解析を進めていくことができます。

しかし「サイト$i,j$の電子って何？」とか、「飛び移り積分って？」とか「電子を作ったり消したりする演算子っていったい何？？？」とか考えてみるとなんだかよくわからない、（学部時代を適当に過ごしてきた）大学院生や、他分野の研究者・技術者の方もいるのではないでしょうか。
というか何を隠そう、私が学部時代に量子力学を適当に勉強した状態で修士に進み、その辺にもやもやを抱えながら研究してきた残念院生だったのですが（そしてそれを解消できないまま修了しちゃったのですが）、このままだと死んでも死にきれないので過去の自分がブラックボックスにしてきた部分を理解するべく勉強した内容をまとめようと思い立ちました。

量子力学についてのなんとなくの理解を前提に、上記Tight-binding modelの第二量子化表示について納得できるまでの内容を順次作成していこうと思います。とはいえ全部書くとすごい量になってしまい終わらなさそうなので、適宜既存の教科書に説明を譲りつつ、教科書では省略されていたりあまり触れられていない部分を中心に書いていこうと思っています。

# 目次（予定）
公開予定の記事は以下の通りです。紙の上では大体できているので、latexを打ったり説明を補足したりでき次第公開予定です。（とはいえ書きながらやっぱりわかっていなかった部分がどんどん出てきて中々進まないのですが。。）

## [水素原子中の電子](https://zenn.dev/ponzumai/articles/tight-binding-model-hydrogen-atom)
## [（補足）水素原子中の電子のシュレディンガー方程式の具体的な計算](https://zenn.dev/ponzumai/articles/appendix-hidrogen-atom)
## [多電子状態の波動関数](https://zenn.dev/ponzumai/articles/tight-binding-model-many-electron)
## [電子のスピンを考慮した多電子系の波動関数](https://zenn.dev/ponzumai/articles/tight-binding-model-spin)
## [多電子原子中の電子](https://zenn.dev/ponzumai/articles/tight-binding-model-many-electron-atom)
## 格子ベクトルと逆格子ベクトル
## 固体（結晶）中の電子状態
## Tight-binding 近似（第一量子化）
## Tight-binding modelの第二量子化表示

# さいごに
質問やわかりにくい点などあればコメントを頂けますと幸いです。（もちろん感想だけでもいただけるととても喜びます。）
また、上述の通り滑り込み修士が独学をベースに書いているノートなので、勘違いや間違いは多数含まれている可能性がありますし、おかしな書き方をしている部分もあるかと思いますのでご注意ください（指摘していただければ嬉しいです）。