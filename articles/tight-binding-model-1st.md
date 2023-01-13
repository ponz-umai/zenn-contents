---
title: "Tight-bindingモデル（第一量子化）"
emoji: "🧵"
type: "tech" # tech: 技術記事 / idea: アイデア
topics: ["quantum","quantumcomputing","quantumcomputer","物理","物理学"]
published: false
---
# はじめに
これまで（水素）原子波動関数から始まり、多電子系の取り扱い、周期ポテンシャルの取り扱いを経て、ようやくTight-bindingモデルを考えるための準備が整いました。

## これまでの振り返り

今まで散々書いてきた内容ですが、改めて私たちが今考えたい固体内の状態から考え始めてみます。

まず、何も近似を入れなければ、固体内の電子は
- 原子核と電子が相互に運動している下で
- 原子核と電子の相互作用があり
- 電子間もそれぞれに相互作用し
- 原子核間も相互作用している

みたいなカオスな状況です。これらの運動エネルギーと相互作用ポテンシャルを全て考慮すると、シュレーディンガー方程式は$3$次元$\times N_e$（電子数）$\times 2$（スピン）$+3$次元$\times N_a$（原子核の数）$\times s_a$（原子核のスピン）個の変数からなる関数の微分方程式となり、全く解けません。

そこで、
- 原子核は静止しており、ついでに内殻電子も原子核の近傍で留まっている→原子核と内殻電子のポテンシャルを全てまとめて、1電子ポテンシャルとして取り入れる
- 価電子間の相互作用も無視する
- スピンに関する相互作用も無視する

ことで、シュレーディンガー方程式は多変数の微分方程式であっても以下のように変数分離が可能な形


$$
\mathcal{H}\Phi(\boldsymbol{r}_1, \boldsymbol{r}_2,\cdots) = E\Phi(\boldsymbol{r}_1, \boldsymbol{r}_2,\cdots),\\
\mathcal{H} = \sum_i\left(  -\frac{\hbar^2}{2m}\nabla_i^2
+
U(\boldsymbol{r}_i) 
\right) 
$$

となります。ここで、色々な相互作用を取り入れた1電子ポテンシャル$U(\boldsymbol{r})$は固体の周期性を表す格子ベクトル

$$
\boldsymbol{n} = n_1\boldsymbol{a}_1 + n_2\boldsymbol{a}_2 + n_3\boldsymbol{a}_3
$$

に対して

$$
U(\boldsymbol{r} + \boldsymbol{n}) = U(\boldsymbol{r})
$$

を満たす周期ポテンシャルで、各電子に対して共通の形を仮定しています。

このように1電子のハミルトニアンの総和で多電子系のハミルトニアンが書かれるとき、[以前確認した](https://zenn.dev/ponzumai/articles/tight-binding-model-spin)ように、1体のハミルトニアンからなるシュレーディンガー方程式

$$
\hat{H}\varphi(\boldsymbol{r}) = \varepsilon\varphi(\boldsymbol{r}),\\

\hat{H} = 
    -\frac{\hbar^2}{2m}\nabla^2
    +
   U(\boldsymbol{r})
$$

を解いて1体固有関数を求め、その固有関数を並べたスレーター行列式（やその線形結合）を考えることで多体の波動関数を得ることができるのでした。

また、[この章で確認したように](https://zenn.dev/ponzumai/articles/tight-binding-model-electrons-in-solids)上記のような周期ポテンシャルの下で、電子が取る固有関数はBloch関数として表され、さらに[こちらで確認したように](https://zenn.dev/ponzumai/articles/tight-binding-model-wannier-func)局在したWannier関数$w_{n,\boldsymbol{R}}(\boldsymbol{r}) = w_n(\boldsymbol{r} - \boldsymbol{R})$の線形結合

$$
\varphi_{n,\boldsymbol{k}}(\boldsymbol{r}) =\sum_{\boldsymbol{R}}w_{n,\boldsymbol{R}}(\boldsymbol{r})e^{i\boldsymbol{k}\cdot\boldsymbol{R}},\\
w_{n,\boldsymbol{R}}(\boldsymbol{r})
=
\frac{1}{N}\sum_{\boldsymbol{k}} \varphi_n(\boldsymbol{r},\boldsymbol{k})e^{-i\boldsymbol{k}\cdot\boldsymbol{R}}
$$

で展開でき、またWannier関数は
- 並進性：$w_{n,\boldsymbol{R}}(\boldsymbol{r}) = w_n(\boldsymbol{r}-\boldsymbol{R})$
- 正規直交性：$\int_V w_{n'}^*(\boldsymbol{r}-\boldsymbol{R}')w_n(\boldsymbol{r} - \boldsymbol{R})dr =\delta_{n,n'}\delta_{\boldsymbol{R},\boldsymbol{R}'}$
- 局在性：ラベル＝格子点の座標$\boldsymbol{R}$を中心として局在した関数

を満たすのでした。

## Tight-bindingモデルの構築

これまでの振り返りを踏まえて、いくつかの仮定や近似を入れながらTight-bindingモデルを構築していきます。

その道のりは、大きく2つの工程に分かれます。

第一に、固体内のポテンシャルを孤立原子ポテンシャルで近似し、かつ固有関数であるBloch関数を原子軌道関数（孤立原子ポテンシャルの固有関数）の線形結合で表す**LCAO法**を用いて、変分法により固体内の電子の1体固有関数を求める方程式を導出します。

第二に、上記方程式やLCAOの展開に対して、「電子が孤立原子に強く束縛されている」と仮定した近似（Tight-binding近似）を行い、具体的な固有関数の形、そして固有エネルギー（エネルギーバンド）を導きます。

**（余談）**この記事を書いていて、色々参考にしていると、同じような操作でも「Tight-binding（又は強束縛、強結合等々）近似」と書いていたり、「Tight-binding（又は以下略）モデル」と書いていたり、どういう部分が「近似」と呼ばれてどこまで行くと「モデル」になるのかよくわからなくなってきました。
その辺の厳密な使い分けがあるのかどうか今のところよくわかっていないのですが、本章では、
- 「電子が強く束縛されている」という描像に基づいた近似を「Tight-binding近似」
- そのような近似によって得られた最終的な波動関数（固有関数）・固有エネルギーを「Tight-bindingモデル」

と呼び分けることにします。全然違ってたら、すみません。


# 第一工程：LCAO＋変分法
さて、「はじめに」で整理した一般的な周期ポテンシャル下、すなわち固体内の電子が満たすべき条件ですが、やはりこれだけではまだまだシュレーディンガー方程式は解けません。

ここから、何とかして解を求められる形に持っていけるように、まずは第一段階として電子の固有関数を具体的な関数形で展開し、その展開係数が満たすべき条件を求めることで、固体内の電子を表すための方程式を獲得します。

## 周期ポテンシャルの近似

すでに[多電子原子中の電子](https://zenn.dev/ponzumai/articles/tight-binding-model-many-electron-atom)の章でちらっと述べましたが、まずは周期ポテンシャル$U(\boldsymbol{r})$を、各格子点の孤立原子（や分子等）＋内殻電子の1体ポテンシャルを並べた関数

$$
U(\boldsymbol{r}) = \sum_{\boldsymbol{R}}V(\boldsymbol{r} - \boldsymbol{R})
$$

で近似します。あるいは、価電子間の相互作用もこのポテンシャルに上手く取り入れることもあるかもしれません、いずれにしても後に求めるエネルギーバンドが、上手く実験結果と合うような1体ポテンシャルを求めて、それが周期的に並んだような固体のポテンシャルで近似します。（この近似に何か特別な名前がついているかどうかは知りません）


具体的には、[多電子原子中の電子](https://zenn.dev/ponzumai/articles/tight-binding-model-many-electron-atom)で用いた多電子原子のHartree Fock近似（HF近似）（≒平均場近似）のように多電子の相互作用を平均したポテンシャルで置き換えたり、あるいは、実際の研究レベルだとさらに進んだ手法が用いられているようです。^[そのあたりは知識不足なので手近な固体物理の教科書を参照してください。「疑ポテンシャル法」とか、密度汎関数理論とか、そういうキーワードをよく目にする気がします。私はまだ怖くてちゃんと勉強していません。]とにかく、何らかの手段で、上手いこと固体内の電子が受けるポテンシャルを格子点に紐づいた関数の総和で書けたとします。

## 孤立原子における電子の固有状態

さらに、上記のように固体内のポテンシャルを上手いこと孤立ポテンシャル（これを総称して「孤立原子ポテンシャル」と呼ぶことにします。考えたい問題によっては孤立「分子」やその他何らかのより適切な表現があるかもしれませんが、ひとまずそういう用語ということで）で置き換えたことで、その孤立原子ポテンシャルに対する固有関数系を求めることが可能となりました。

具体的には、例えば[多電子原子中の電子](https://zenn.dev/ponzumai/articles/tight-binding-model-many-electron-atom)の章の後半で述べたように、HF近似により水素様原子のシュレディンガー方程式のポテンシャル部分に、球対称な平均場ポテンシャルが付け加わった形


$$
V(\boldsymbol{r}) = -\frac{Ze^2}{4\pi\epsilon _0r}  
 + \sum_lV_l^c(r) + \sum_lV_l^{ex}(r)
$$

を**局所1体ポテンシャル**として考えたとします。
ここで$V_l^c$は$l$番目の電子に関するCoulomb積分、$V_l^{ex}$は$l$番目の電子との交換積分（を球対称化したもの）です。詳細は上記リンクをご確認ください（しつこいようですがもっと進んだ1体ポテンシャルを使う場合には適宜上式のポテンシャル部分を書き換えたようなものを考えることになるでしょう）。

このようなポテンシャルの下での孤立原子シュレーディンガー方程式

$$
\left(-\frac{\hbar^2}{2m} \nabla^2  -\frac{Ze^2}{4\pi\epsilon _0r}  
 + \sum_lV_l^c(r) + \sum_lV_l^{ex}(r)

\right) 

\phi_m(\boldsymbol{r}) =\varepsilon_m^{\rm a} \phi_m (\boldsymbol{r} )
$$

の下で求めた固有エネルギーを「原子準位」と呼んだりします。また以降、atomの"a"を取って、$m$番目の固有エネルギーを$\varepsilon_m^{\rm a}$と表し（エネルギーの小さい順からラベル$n$が振られているとします）、その固有エネルギーに対応した固有関数を$\phi_m(\boldsymbol{r})$とします。

もう少し孤立原子の固有状態について考えていくと、前掲の章で述べたように、方程式は球対称ポテンシャルの近似のもと以下のように変数分離ができ、

$$
-\frac{\hbar^2}{2m} \left(
        \frac{d ^2R}{d r^2} + \frac{2}{r}\frac{d R}{d r}
        -\frac{\lambda }{r^2} R
        \right)
        +
        \left(
        -\frac{Ze^2}{4\pi\epsilon _0r}
        + \sum_lV_l^c(r) + \sum_lV_l^{ex}(r)
        \right)R = \varepsilon R
       ,\\
\frac{1}{\sin\theta } \frac{d}{d\theta } \left( \sin\theta \frac{d\Theta }{d\theta }   \right) + \left( \lambda - \frac{m^2}{\sin^2\theta }  \right) \Theta =0,\\
\frac{d^2\Phi }{d\phi ^2} + m^2\Phi =0.
$$

角度部分は水素様原子とまったく同じ形で固有関数は球面調和関数

$$
Y_l^m(\theta,\phi)
$$

で書くことができます。

ただ動径関数部分は水素原子の固有関数とは異なるため水素原子の固有状態とは異なるものとはなりますが、動径部分の固有状態を何かしら求めたとして、その固有値（固有エネルギー）の大きさに応じてラベル$n$でラベル付けし、それと球面調和関数の積を固有状態と考えることで**水素原子の固有状態で用いた$1s,2s,2p_x,\cdots$状態と類似したもの**と考えることができます。


## 原子軌道関数の線形結合を用いたWannier関数の近似（LCAO近似） 

上記のようにして求めた固有関数形$\{\phi_m\}$は（多分）^[孤立原子シュレーディンガー方程式の固有関数系が完全系となることについてちょっと自信がないのですが、例えば動径関数を何らかの完全系$\{\varphi_n(r)\}$の線形結合で、角度部分は球面調和関数で展開したとして、そのように展開された固有関数系を全て（無限個）集めたとしたら、それは結局行列表示した孤立原子ハミルトニアンを対角化する行列による、完全系$\{\psi_l(\boldsymbol{r})\}( = \{\varphi_n(r) Y_l^m(\theta,\phi)\})$のユニタリ変換とみなせて、それは結局完全系を満たす、というようなことなのかと思っています。なんか出鱈目なことを書いていそうで怖いですが]完全系を成すので、**固体中の電子の波動関数を原子軌道関数で展開**することができます。

特に、前章で求めたように、固体内の電子は「格子点に局在した関数」Wannier関数の線形結合で表現できるのでした。このWannier関数を、孤立原子（又は分子等）の波動関数$\phi_m(\boldsymbol{r})$と展開係数$b_m^n$を用いて、

$$
w_{n,\boldsymbol{R}}(\boldsymbol{r}) = w_n(\boldsymbol{r}-\boldsymbol{R}) \simeq \sum_mb_m^n\phi_m(\boldsymbol{r}-\boldsymbol{R}) 
$$

と展開します。（原子軌道関数であることがわかりやすいように、以降、原子軌道関数は先ほどと同じく$\phi$を用いて表します。）

このような考え方は、「原子軌道の線形結合」の英語バージョン"Linear Combination of Atomic Orbital"の頭文字をとってLCAO法、またはLCAO近似（展開する原子軌道に制限をつけなければ「近似」ではないですが、多くの場合物理的な考察に応じて少数の原子軌道で展開することになります。（そうしないとそもそも計算できないし）そのように少数の原子軌道で展開するときに「LCAO近似」というのだと思います）と呼ばれています。

この結果、結晶内の電子の固有状態＝Bloch関数は

$$
\begin{align*}
\varphi_{n,\boldsymbol{k}} &= \sum_{\boldsymbol{R}}e^{i\boldsymbol{k}\cdot\boldsymbol{R}}w_n(\boldsymbol{r}-\boldsymbol{R})\\

&=
 \sum_{\boldsymbol{R}}e^{i\boldsymbol{k}\cdot\boldsymbol{R}} \sum_mb_m^n\phi_m(\boldsymbol{r}-\boldsymbol{R})
 

\end{align*}
$$

と書くことができます。（一応、展開に使う関数にまだ制限は掛けていないので、まだ等式で書いています）

このように固体内の固有関数を、既知の関数系$\{\phi_m\}$で表現することができました。この展開を出発点として、展開係数を求めることができれば固有状態・固有エネルギーが分かります。

::: details （余談その2）
余談ですが、LCAO近似とTight-binding近似は割と混同されがちですが、別物であるようで、本章でもそれを意識した書き方にしています。
これは多分、原子軌道の線形結合を取ったとしても、十分空間的に広がった（束縛が弱い）電子の波動関数を表すこともできるし、また逆に強く束縛された電子の関数を表現するための方法は必ずしも原子軌道の線形結合でなくても良い、という事なのだと思いますが、詳細はわかりません。
とはいえ、ここで書いたことはこれから述べていく近似を追っていくことで具体的にイメージできるかと思います。
:::


## 変分法を用いた固有値方程式の導出

### 変分法の超概要

さて、Wannier関数を原子軌道を用いて展開しました。

$$
\varphi_{n,\boldsymbol{k}}

=
 \sum_{\boldsymbol{R}}e^{i\boldsymbol{k}\cdot\boldsymbol{R}} \sum_mb_m^n\phi_m(\boldsymbol{r}-\boldsymbol{R})
 
$$

ここから、上記関数を試行関数として、具体的な関数形、つまり係数$b_m^n$を求めることができれば、具体的な固有関数、そして固有エネルギーを求めることができます。

この係数を決定するために、過去にも用いた「変分法」を使います。すなわち、エネルギー期待値を最小化するような展開係数を求めていきます。

一旦詳細を省略しますが、「（試行）関数$\varphi_{n,\boldsymbol{k}}(\boldsymbol{r})$がエネルギー期待値を最小にする」という条件式は結局、元々考えていた固体内の1体ハミルトニアンに対するシュレーディンガー方程式

$$
\hat{H}\varphi_{n,\boldsymbol{k}}(\boldsymbol{r}) = \varepsilon_{n,\boldsymbol{k}}\varphi_{n,\boldsymbol{k}}(\boldsymbol{r}),\\

\hat{H} = 
    -\frac{\hbar^2}{2m}\nabla^2
    +
   \sum_{\boldsymbol{R}}V(\boldsymbol{r} - \boldsymbol{R})
$$

と一致します。ここでハミルトニアンが孤立原子のハミルトニアンではなく、孤立原子のポテンシャルを周期的に並べた固体のハミルトニアンであることに注意してください（念のため）

上式に$\varphi_{n,\boldsymbol{k}}=\sum_{\boldsymbol{R}}e^{i\boldsymbol{k}\cdot\boldsymbol{R}} \sum_mb_m^n\phi_m(\boldsymbol{r}-\boldsymbol{R})$を代入し、

$$
\left(
 -\frac{\hbar^2}{2m}\nabla^2
    +
   \sum_{\boldsymbol{R}}V(\boldsymbol{r} - \boldsymbol{R})
   \right)
   \sum_{\boldsymbol{R}}e^{i\boldsymbol{k}\cdot\boldsymbol{R}} \sum_mb_m^n\phi_m(\boldsymbol{r}-\boldsymbol{R})
    = \varepsilon_{n,\boldsymbol{k}}\sum_{\boldsymbol{R}}e^{i\boldsymbol{k}\cdot\boldsymbol{R}} \sum_mb_m^n\phi_m(\boldsymbol{r}-\boldsymbol{R})
$$

原点を中心とした原子軌道関数の複素共役$\phi_{l}^*(\boldsymbol{r})$（$\boldsymbol{R} = \boldsymbol{0}$の場合に対応する）をかけて全空間について積分することで係数を決定するための固有値方程式を導いていきます。




### 固有値方程式の導出



・・・とはいえ、これでは近似も何もなくただのシュレーディンガー方程式で、このままで解けるなら何の苦労もないわけです。そこでまず初めに、重要な関係性として**Wannier関数は少数の原子軌道関数で精度良く展開できる**ことを確かめます。^[私はこの部分が最近までイマイチ理解できておらず、Tight-bindingモデルの理解にかなり遠回りをしました]


### 展開に用いる原子軌道の見積もり^[この部分はアシュクラフト・マーミン上（I)を参考にした。]^[ものすごく余談ですが、最初にゼミで読んだ翻訳書がニジェ・オーランドでGrassmann代数をさらっとかじった後アルトランド・サイモンズで経路積分をかじってさらに例題を解いてくみたいな内容だったので、カナカナ連名著者の教科書は自分には敷居の高い高尚な内容が書かれているものだと思い込んでアシュクラフト・マーミン等も避けていたのですが、最近になって読んでみると非常に丁寧に知りたいことが書いてあり愕然としました。さっさと読めばよかった。。。]

さて、Wannier関数を原子軌道関数で展開したのですが、展開に用いる原子軌道が多くなればなるほど精度が良くなるように思えます。それはそうなのですが、そんなことを言うと近似も近似でなくなるわけで、以下のようにある係数$b_m$、つまり展開に含まれる原子軌道のうち、$m$番目のエネルギー順位の原子軌道が含まれる割合についての関係式を得ることができ、結局そのエネルギーバンド程度の固有エネルギーに対応する（少数の）原子軌道関数のみでWannier関数を展開できることが分かります。

早速前に進みましょう。先ほど変分法によって得られた、固有状態が満たすべき方程式（＝シュレーディンガー方程式）

$$
\left(
 -\frac{\hbar^2}{2m}\nabla^2
    +
   \sum_{\boldsymbol{R}}V(\boldsymbol{r} - \boldsymbol{R})
   \right)
   \varphi_{n,\boldsymbol{k}}(\boldsymbol{r}) = \varepsilon_{n,\boldsymbol{k}}\varphi_{n,\boldsymbol{k}}(\boldsymbol{r})  \\
   
\Rightarrow
\left(
 -\frac{\hbar^2}{2m}\nabla^2
    +
   \sum_{\boldsymbol{R}}V(\boldsymbol{r} - \boldsymbol{R})
   \right)
   \sum_{\boldsymbol{R}}e^{i\boldsymbol{k}\cdot\boldsymbol{R}} \sum_mb_m^n\phi_m(\boldsymbol{r}-\boldsymbol{R})
    = \varepsilon_{n,\boldsymbol{k}}\sum_{\boldsymbol{R}}e^{i\boldsymbol{k}\cdot\boldsymbol{R}} \sum_mb_m^n\phi_m(\boldsymbol{r}-\boldsymbol{R})
$$

に、原点を中心とした原子軌道関数の複素共役$\phi_{l}^*(\boldsymbol{r})$（$\boldsymbol{R} = \boldsymbol{0}$の場合に対応する）をかけて全空間について積分します。

まずは比較的シンプルな右辺について、右辺後の計算がわかりやすいように$\sum_{\boldsymbol{R}}$を$\boldsymbol{R} = \boldsymbol{0}$とそれ以外に分けて書くと、

$$
\begin{align*}
（右辺の積分）&=\varepsilon_{n,\boldsymbol{k}}
\left(
   \sum_m b^n_m \int \phi_{l}^*(\boldsymbol{r})\phi_{m}(\boldsymbol{r})d\boldsymbol{r}

   +
   \sum_{\boldsymbol{R}\neq\boldsymbol{0}}\sum_m b^n_m e^{i\boldsymbol{k}\cdot\boldsymbol{R}}
   
   \int \phi_{l}^*(\boldsymbol{r})\phi_{m}(\boldsymbol{r} - \boldsymbol{R})d\boldsymbol{r}

\right)\\

&=
\varepsilon_{n,\boldsymbol{k}}
\left(
    b^n_l
   +
   \sum_m 
   \left[
   \sum_{\boldsymbol{R}\neq\boldsymbol{0}}
   e^{i\boldsymbol{k}\cdot\boldsymbol{R}}
   
   \int \phi_{l}^*(\boldsymbol{r})\phi_{m}(\boldsymbol{r} - \boldsymbol{R})d\boldsymbol{r}
   \right]b^n_m 

\right)

\end{align*}
$$

となります。ここで2式目は、原子軌道関数の正規直交性

$$
\int \phi_{l}^*(\boldsymbol{r})\phi_{m}(\boldsymbol{r})d\boldsymbol{r} = \delta_{l,m}
$$

を使いました。

また、左辺については、原子軌道関数が固体中の1体ハミルトニアンのうち、原点以外の局所ポテンシャルを除いた部分の固有関数であること：

$$
\left(-\frac{\hbar^2}{2m} \nabla^2  + V(\boldsymbol{r})

\right) 

\phi_l(\boldsymbol{r}) =\varepsilon_l^{\rm a} \phi (\boldsymbol{r} )
$$

を思い出し、またエルミート演算子であるハミルトニアンの積分に関する性質：

$$
\begin{align*}
\int \varphi^*(\boldsymbol{r})
\left (
   \hat{H} \psi(\boldsymbol{r}) 
\right) d\boldsymbol{r}
&= \int \left(\hat{H}\varphi(\boldsymbol{r})\right)^* \psi(\boldsymbol{r})
\end{align*}
$$

特に$\varphi$が$\hat{H}$の固有状態のとき、固有値を$\varepsilon$とすると固有値は実数なので、$\int \left(\hat{H}\varphi(\boldsymbol{r})\right)^* \psi(\boldsymbol{r}) = \varepsilon\int \varphi^*(\boldsymbol{r}) \psi(\boldsymbol{r})$、を利用して積分を計算していきます。

こちらも$\boldsymbol{R} = \boldsymbol{0}$の場合とそれ以外、さらにハミルトニアン部分も孤立原子のハミルトニアンとそれ以外で分けて

$$
\begin{align*}
(左辺の積分)&=
\sum_{\boldsymbol{R}}\sum_me^{i\boldsymbol{k}\cdot\boldsymbol{R}} b_m^n
\int \phi_l^*(\boldsymbol{r})
\left(
 -\frac{\hbar^2}{2m}\nabla^2
    +
   \sum_{\boldsymbol{R}'}V(\boldsymbol{r} - \boldsymbol{R}')
   \right)
   \phi_m(\boldsymbol{r}-\boldsymbol{R}) d\boldsymbol{r}\\

&=
\sum_m \left[
    \int \phi_l^*(\boldsymbol{r})
   \left(
   -\frac{\hbar^2}{2m}\nabla^2
      +
      V(\boldsymbol{r})
      \right) \phi_m(\boldsymbol{r})d\boldsymbol{r}
      \right] b^n_m \\

&\>\>\>\>+
\sum_m
\left[ 
   \sum_{\boldsymbol{R}\neq\boldsymbol{0}}
   e^{i\boldsymbol{k}\cdot\boldsymbol{R}} 
  \int \phi_l^*(\boldsymbol{r})
   \left(
   -\frac{\hbar^2}{2m}\nabla^2
      +
      V(\boldsymbol{r} )
      \right) \phi_m(\boldsymbol{r}-\boldsymbol{R})d\boldsymbol{r}
      \right]  b^n_m 
      \\

&\>\>\>\>+
\sum_m
\left[
    \int \phi_l^*(\boldsymbol{r})
   \left(
   \sum_{\boldsymbol{R}'\neq \boldsymbol{0}}
      V(\boldsymbol{r}-\boldsymbol{R}')
      \right) \phi_m(\boldsymbol{r})d\boldsymbol{r}
      \right] b^n_m 
      \\


&\>\>\>\>+
\sum_m 
\left[
   \sum_{\boldsymbol{R}\neq\boldsymbol{0}}
   e^{i\boldsymbol{k}\cdot\boldsymbol{R}}  \int \phi_l^*(\boldsymbol{r})
   \left(
   \sum_{\boldsymbol{R}'\neq \boldsymbol{0}}
      V(\boldsymbol{r}-\boldsymbol{R}')
      \right)  \phi_m(\boldsymbol{r}-\boldsymbol{R})d\boldsymbol{r}
      \right] b^n_m\\

&=
\varepsilon_l^{\rm a}  b^n_l \\

&\>\>\>\>+
\varepsilon_l^{\rm a}
\sum_m 
\left[
   \sum_{\boldsymbol{R}\neq\boldsymbol{0}}
   e^{i\boldsymbol{k}\cdot\boldsymbol{R}} 
    \int \phi_l^*(\boldsymbol{r})
   \phi_m(\boldsymbol{r}-\boldsymbol{R})d\boldsymbol{r}
   \right] b^n_m\\



&\>\>\>\>+
\sum_m 
\left[
    \int \phi_l^*(\boldsymbol{r})
   \left(
   \sum_{\boldsymbol{R}'\neq \boldsymbol{0}}
      V(\boldsymbol{r}-\boldsymbol{R}')
      \right) \phi_m(\boldsymbol{r})d\boldsymbol{r}
      \right] b^n_m\\

&\>\>\>\>+
\sum_m 
\left[
   \sum_{\boldsymbol{R}\neq\boldsymbol{0}}
   e^{i\boldsymbol{k}\cdot\boldsymbol{R}} \int \phi_l^*(\boldsymbol{r})
   \left(
   \sum_{\boldsymbol{R}'\neq \boldsymbol{0}}
      V(\boldsymbol{r}-\boldsymbol{R}')
      \right)  \phi_m(\boldsymbol{r}-\boldsymbol{R})d\boldsymbol{r}
      \right]b^n_m  

\end{align*}
$$

となります。3式目の第1項、第2項の変形に、先ほど述べた原子軌道関数が孤立原子の固有状態であること、および固有値が実数であることを使いました。

さて、先ほどの等式から$b^n_l$単体の項とそれ以外をそれぞれ移項したりして整理すると、展開係数$b^n_l$に関する等式

$$
\begin{align*}
(\varepsilon_{n,\boldsymbol{k}} - \varepsilon_l^{\rm a})b_l^n 

&=
-(\varepsilon_{n,\boldsymbol{k}} - \varepsilon_l^{\rm a})
\sum_m 
\left[
   \sum_{\boldsymbol{R}\neq\boldsymbol{0}}
   e^{i\boldsymbol{k}\cdot\boldsymbol{R}} 
   \int \phi_l^*(\boldsymbol{r})
   \phi_m(\boldsymbol{r}-\boldsymbol{R})d\boldsymbol{r}
   \right]b^n_m \\

&\>\>\>\>+
\sum_m 
\left[
    \int \phi_l^*(\boldsymbol{r})
   \left(
   \sum_{\boldsymbol{R}'\neq \boldsymbol{0}}
      V(\boldsymbol{r}-\boldsymbol{R}')
      \right) \phi_m(\boldsymbol{r})d\boldsymbol{r}
      \right] b^n_m\\

&\>\>\>\>+
\sum_m 
\left[
   \sum_{\boldsymbol{R}\neq\boldsymbol{0}}
   e^{i\boldsymbol{k}\cdot\boldsymbol{R}} \int \phi_l^*(\boldsymbol{r})
   \left(
   \sum_{\boldsymbol{R}'\neq \boldsymbol{0}}
      V(\boldsymbol{r}-\boldsymbol{R}')
      \right)  \phi_m(\boldsymbol{r}-\boldsymbol{R})d\boldsymbol{r}
      \right]b^n_m 
\end{align*}
$$

を得ることができます。

等式の右辺についてひとつづつ見ていくと、まず第1項は異なる場所に中心を持つ原子軌道関数同士の積分が含まれています。ここで今、「格子点に強く束縛された状態」を想定していたことから、その展開に含まれる原子軌道関数も各格子点に局在していると考え、この積分の値は小さいとみなせます。

また第2項についても、原点を中心とする原子軌道関数$\phi(\boldsymbol{r})$は格子間距離$\boldsymbol{R}$離れた場所では小さくなる一方、原点においては原点以外を中心として持つ孤立原子ポテンシャル$V(\boldsymbol{r}-\boldsymbol{R})$が小さくなると考えられるので、こちらの積分も小さい値を持つと考えられます。
^[アシュクラフト・マーミンでは「孤立原子ポテンシャルは波動関数と比べるとそこまで急速に小さくなるわけではないのでこの項が「小さい」とするのは幾分怪しいが、$\boldsymbol{k}$に依存しない項なので問題ないのと、ポテンシャルと波動関数を適切に選べばこの項も小さくできる」というようなことが脚注で書いてある。$\boldsymbol{k}$に依存しないとなぜ問題ないのか、どうやって小さくするのか等まだよくわかっていないのだがひとまず受け入れておくことにする]

最後に第3項についても、第1項、第2項と似たような理由で、中心が異なる（原点、$\boldsymbol{R}'$、$\boldsymbol{R}$）3つの局在した関数の積分となるので値は小さいと考えられます（これを「3中心積分などと呼んだりする」）。

以上のように、右辺はどの項も小さい値を持つと見積もることができ、左辺

$$
(\varepsilon_{n,\boldsymbol{k}} - \varepsilon_l^{\rm a})b_l^n 
$$

は小さくならなければならないことがわかります。これを満たすには

$$
\varepsilon_{n,\boldsymbol{k}} \simeq \varepsilon_l^{\rm a}
$$

すなわち今考えているバンドのエネルギーに対して、係数$b_l^n$に対応した原子軌道関数の固有エネルギー$\varepsilon_l^{\rm a}$が近い値を持つ場合、または

$$
b_l^n \simeq 0
$$

の場合であると言えます。

[固体（結晶）中の電子状態](https://zenn.dev/ponzumai/articles/tight-binding-model-electrons-in-solids)の章で見たように、固体内の電子の固有状態であるBloch関数$\varphi_{n,\boldsymbol{k}}$は、
- バンド指標$n$ごとに離散的な、すなわち程度離れたエネルギーを持ち
- エネルギー$\varepsilon_{n,\boldsymbol{k}}$は連続、かつ限定的な領域「ブリルアンゾーン」を周期として持つ周期関数

なのでした。つまり指標$n$の固有エネルギーは




従って、「あるバンド$n$のBloch関数$\varphi_{n,\boldsymbol{k}}(\boldsymbol{r})$の展開に用いるWannier関数$w_{n,\boldsymbol{R}}(\boldsymbol{r})$の展開に用いる原子軌道関数$\phi_l(\boldsymbol{r})$」は、そのBloch関数が属するエネルギーバンド


### Tight-bindingな近似その1

さて、上記LCAO近似で電子の固有状態を近似しましたが、ここでより「強く束縛された状態」の描像に基づいた近似を行います。

具体的には、Wannier関数の近似に用いた原子軌道$varphi_m(\boldsymbol{r}) $に対して、異なる格子点にある原子軌道同士の「重なり」が小さいとして、次の積分を$0$と考えます。つまり


$$
\int\varphi_m^*(\boldsymbol{r}) \varphi_l(\boldsymbol{r-\boldsymbol{R}})d\boldsymbol{r} \simeq 0 
$$

とします。この積分は **Overlap Integral（重なり積分）** と呼ばれ、格子間距離に比べてWannier関数を近似している原子軌道の広がりが十分小さい、つまり「強く束縛されている」ようなイメージに対応しています。

また、これによりLCAO近似で表したWannier関数も正規直交性を持ちます。（Bloch関数のFourier展開として定義されるWannier関数は正規直交性を満たしますが、


### Tight-bindingな近似その2

隣接格子

### 固有状態と固有エネルギー








また先ほど記号を定義した孤立原子の固有エネルギーを、対応する原子軌道関数の状態$s,p,d,\cdots$に応じて$\varepsilon_s^{\rm a}, \varepsilon_p^{\rm a}, \cdots$等と書くこともあります。この時「主量子数」$n$についてはあまり言及されないことが多いですが、これは後に出てくるように、Tight-binding近似を行うにあたりその角度部分の性質（符号等）がバンドの形状に大きな影響を与え、軌道部分は