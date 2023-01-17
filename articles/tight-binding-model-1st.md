---
title: "第一量子化のTight-bindingモデル"
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

また、[この章で確認したように](https://zenn.dev/ponzumai/articles/tight-binding-model-electrons-in-solids)上記のような周期ポテンシャルの下で、電子が取る固有関数はBloch関数として表され、さらに[こちらで確認したように](https://zenn.dev/ponzumai/articles/tight-binding-model-wannier-func)局在したWannier関数$w_{\boldsymbol{R}}(\boldsymbol{r}) = w(\boldsymbol{r} - \boldsymbol{R})$の線形結合

$$
\varphi_{\boldsymbol{k}}(\boldsymbol{r}) =\sum_{\boldsymbol{R}}w_{\boldsymbol{R}}(\boldsymbol{r})e^{i\boldsymbol{k}\cdot\boldsymbol{R}},\\
w_{\boldsymbol{R}}(\boldsymbol{r})
=
\frac{1}{N}\sum_{\boldsymbol{k}} \varphi_{\boldsymbol{k}}(\boldsymbol{r})e^{-i\boldsymbol{k}\cdot\boldsymbol{R}}
$$

で展開でき、またWannier関数は
- 並進性：$w_{\boldsymbol{R}}(\boldsymbol{r}) = w(\boldsymbol{r}-\boldsymbol{R})$
- 正規直交性：$\int_V w_{n'}^*(\boldsymbol{r}-\boldsymbol{R}')w_n(\boldsymbol{r} - \boldsymbol{R})dr =\delta_{n,n'}\delta_{\boldsymbol{R},\boldsymbol{R}'}$
- 局在性：ラベル＝格子点の座標$\boldsymbol{R}$を中心として局在した関数

を満たすのでした。

## Tight-bindingモデルの構築について

これまでの振り返りを踏まえて、いくつかの仮定や近似を入れながらTight-bindingモデルを構築していきます。

その道のりは、大きく2つの工程に分かれます。

第一に、固体内のポテンシャルを孤立原子ポテンシャルで近似し、かつ固有関数であるBloch関数を原子軌道関数（孤立原子ポテンシャルの固有関数）の線形結合で表す**LCAO法（近似）**を用いて、変分法により固体内の電子の1体固有関数を求める方程式を導出します。

第二に、上記方程式やLCAOの展開に対して、「電子が孤立原子に強く束縛されている」と仮定した近似（Tight-binding近似）を行い、具体的な固有関数の形、そして固有エネルギー（エネルギーバンド）を導きます。

**（余談）** この記事を書いていて、色々参考にしていると、同じような操作でも「Tight-binding（又は強束縛、強結合等々）近似」と書いていたり、「Tight-binding（又は以下略）モデル」と書いていたり、どういう部分が「近似」と呼ばれてどこまで行くと「モデル」になるのかよくわからなくなってきました。
その辺の厳密な使い分けがあるのかどうか今のところよくわかっていないのですが、本章では、
- 「電子が強く束縛されている」という描像に基づいた近似を「Tight-binding近似」
- そのような近似によって得られた最終的な固有エネルギー（エネルギーバンド）を「Tight-bindingモデル」

と呼び分けることにします。全然違ってたら、すみません。


# 第一工程：周期ポテンシャルの近似＋LCAO法＋変分法
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
w_{\boldsymbol{R}}(\boldsymbol{r}) = w(\boldsymbol{r}-\boldsymbol{R}) = \sum_mb_m\phi_m(\boldsymbol{r}-\boldsymbol{R}) 
$$

と展開します。（原子軌道関数であることがわかりやすいように、以降、原子軌道関数は先ほどと同じく$\phi$を用いて表します。）ここで$b_m^n$は規格化条件$\sum_m|b_m|^2 = 1$を満たすものとします。

特に、今回のように強く束縛されている系であれば、**固体内の電子状態は概ね孤立原子の束縛状態に近いはず**と考え、有限個の原子軌道での展開を考えます。


このような考え方は、「原子軌道の線形結合」の英語バージョン"Linear Combination of Atomic Orbitals"の頭文字をとってLCAO法、またはLCAO近似（展開する原子軌道に制限をつけなければ「近似」ではないですが、多くの場合物理的な考察に応じて少数の原子軌道で展開することになります。（そうしないとそもそも計算できないし）そのように少数の原子軌道で展開するときに「LCAO近似」というのだと思います）と呼ばれています。

この結果、結晶内の電子の固有状態＝Bloch関数は

$$
\begin{align*}
\varphi_{\boldsymbol{k}} &= \sum_{\boldsymbol{R}}e^{i\boldsymbol{k}\cdot\boldsymbol{R}}w(\boldsymbol{r}-\boldsymbol{R})\\

&\simeq
 \sum_{\boldsymbol{R}}e^{i\boldsymbol{k}\cdot\boldsymbol{R}} \sum_mb_m\phi_m(\boldsymbol{r}-\boldsymbol{R})
 

\end{align*}
$$

と近似されます。特に$\sum_m$は有限個の和を意味します。

このように固体内の固有関数を、既知の関数$\{\phi_m\}$で表現することができました。この展開を出発点として、展開係数を求めることができれば固有状態・固有エネルギーが分かります。

::: details （余談その2）
余談ですが、LCAO近似とTight-binding近似は割と混同されがちですが、別物であるようで、本章でもそれを意識した書き方にしています。
これは多分、原子軌道の線形結合を取ったとしても、十分空間的に広がった（束縛が弱い）電子の波動関数を表すこともできるし、また逆に強く束縛された電子の関数を表現するための方法は必ずしも原子軌道の線形結合でなくても良い、という事なのだと思いますが、詳細はわかりません。
とはいえ、ここで書いたことはこれから述べていく近似を追っていくことで具体的にイメージできるかと思います。
:::


## 変分法を用いた固有値方程式の導出

### 変分法の超概要

さて、Bloch関数を原子軌道を用いて以下のように近似できると考えました。

$$
\varphi_{\boldsymbol{k}}(\boldsymbol{r})

\simeq
 \sum_{\boldsymbol{R}}e^{i\boldsymbol{k}\cdot\boldsymbol{R}} \sum_mb_m\phi_m(\boldsymbol{r}-\boldsymbol{R})
 
$$

あとは展開係数を求めれば固有関数を近似的に求めることができるわけです。どうやって求めるかというと、多電子原子の章で[He原子の固有状態を求めた](https://zenn.dev/ponzumai/articles/tight-binding-model-many-electron-atom)際に用いた変分法の考え方を用います。

一旦証明は省略しますが、He原子の際は試行関数として2電子のスレーター行列式を用いたことから変分法の結果平均場を含む方程式が出てきましたが、今回のような1電子のハミルトニアンの下では、「（試行）関数$\varphi_{n,\boldsymbol{k}}(\boldsymbol{r})$がエネルギー期待値を最小にする」という条件式は結局、元々考えていた固体内の1体ハミルトニアンに対するシュレーディンガー方程式

$$
\hat{H}\varphi_{\boldsymbol{k}}(\boldsymbol{r}) = \varepsilon_{\boldsymbol{k}}\varphi_{\boldsymbol{k}}(\boldsymbol{r}),\\

\hat{H} = 
    -\frac{\hbar^2}{2m}\nabla^2
    +
   \sum_{\boldsymbol{R}}V(\boldsymbol{r} - \boldsymbol{R}),\\

\varphi_{\boldsymbol{k}}=\sum_{\boldsymbol{R}}e^{i\boldsymbol{k}\cdot\boldsymbol{R}} \sum_mb_m\phi_m(\boldsymbol{r}-\boldsymbol{R})
$$

と一致します。ここでハミルトニアンが孤立原子のハミルトニアンではなく、孤立原子のポテンシャルを周期的に並べた固体のハミルトニアンであることに注意してください（念のため）

式を全部具体的に書くとこんな感じになります：

$$
\left(
 -\frac{\hbar^2}{2m}\nabla^2
    +
   \sum_{\boldsymbol{R}}V(\boldsymbol{r} - \boldsymbol{R})
   \right)
   \sum_{\boldsymbol{R}}e^{i\boldsymbol{k}\cdot\boldsymbol{R}} \sum_mb_m\phi_m(\boldsymbol{r}-\boldsymbol{R})
    = \varepsilon_{\boldsymbol{k}}\sum_{\boldsymbol{R}}e^{i\boldsymbol{k}\cdot\boldsymbol{R}} \sum_mb_m\phi_m(\boldsymbol{r}-\boldsymbol{R})
$$



というわけでこの「シュレーディンガー方程式」を満たすBloch関数を求めれば良いわけです。特にこの方程式から離散的な$n$でラベルされた展開に用いた個数の固有値$\varepsilon_{n,\boldsymbol{k}}$と、その値に対応した固有関数

$$
\varphi_{n\boldsymbol{k}}(\boldsymbol{r})

=
 \sum_{\boldsymbol{R}}e^{i\boldsymbol{k}\cdot\boldsymbol{R}} \sum_mb_m^n\phi_m(\boldsymbol{r}-\boldsymbol{R})
 
$$

が得られます。なお、変分法の操作は「エネルギー期待値を最小にする」固有関数と、その時の（近似的な）固有エネルギーを求める、というものでした。そうでありながら、上記方程式を解いて得られた、最低値以外の固有値・固有状態も、通常のシュレーディンガー方程式と同様に2番目、3番目の固有値・固有状態とみなしていいのか？という疑問がありますが、近似的に得られた$i$番目の固有値を$\varepsilon_i^v$（変分"variation"の"v"）、真の固有値を$\varepsilon_i$として、それぞれ$\varepsilon_2 \leq \varepsilon_2^v, \varepsilon_3 \leq \varepsilon_3^v, \cdots \varepsilon_i \leq \varepsilon_i^v, \cdots$が成り立つようです（これも証明は一旦省きます）ので、安心して近似的に得られた固有値をそのまま各バンドの固有エネルギーとして良さそうです。


さて、具体的な固有値方程式を求めるにあたり、ここで、Bloch関数の展開$\varphi_{\boldsymbol{k}}=\sum_{\boldsymbol{R}}e^{i\boldsymbol{k}\cdot\boldsymbol{R}} \sum_mb_m\phi_m(\boldsymbol{r}-\boldsymbol{R})$で和の順序を入れ替えてみると、

$$
\varphi_{\boldsymbol{k}}(\boldsymbol{r})=
\sum_mb_m
\left[\sum_{\boldsymbol{R}}e^{i\boldsymbol{k}\cdot\boldsymbol{R}} \phi_m(\boldsymbol{r}-\boldsymbol{R})
\right]
$$

となり、これはBloch関数を展開係数$b_m$、関数$\sum_{\boldsymbol{R}}e^{i\boldsymbol{k}\cdot\boldsymbol{R}} \phi_m(\boldsymbol{r}-\boldsymbol{R})$の線形結合で展開していると見ることもできます。
この関数はBloch sum（ブロッホ和、Bloch和）と呼ばれるみたいで、Blochの定理を満たします。また、一般に異なる格子点を中心とした原子軌道関数は直交しないので、異なる原子軌道からなるブロッホ和も直交しません。^[正規直交性を満たすようなブロッホ和の構成方法もあるみたいなのですが、この章では必要にならないので今はスルーしておきます。]

式を簡略化するためにBloch和を

$$
\begin{align*}
\varphi_{\boldsymbol{k}}&=
\sum_mb_m
\left[\sum_{\boldsymbol{R}}e^{i\boldsymbol{k}\cdot\boldsymbol{R}} \phi_m(\boldsymbol{r}-\boldsymbol{R})
\right]

&\equiv
\sum_mb_m\Phi_{m,\boldsymbol{k}}(\boldsymbol{r})

\end{align*}
$$

と置き、上記「シュレーディンガー方程式」に左からBloch和$\Phi_n^*$を掛けて全空間について積分し、

$$
\sum_mb_m\int \Phi_{n,\boldsymbol{k}}^*(\boldsymbol{r})\hat{H}\Phi_{m,\boldsymbol{k}}(\boldsymbol{r})d\boldsymbol{r}
=
\varepsilon_{\boldsymbol{k}}\sum_mb_m\int \Phi_{n,\boldsymbol{k}}^*(\boldsymbol{r})\Phi_{m,\boldsymbol{k}}(\boldsymbol{r})d\boldsymbol{r}
$$

行列

$$
\left(M_{\boldsymbol{k}}\right)_{nm} \equiv \int \Phi_{n,\boldsymbol{k}}^*(\boldsymbol{r})\hat{H}\Phi_{m,\boldsymbol{k}}(\boldsymbol{r})d\boldsymbol{r},\\
\left(S_{\boldsymbol{k}}\right)_{nm} \equiv \int \Phi_{n,\boldsymbol{k}}^*(\boldsymbol{r})\Phi_{m,\boldsymbol{k}}(\boldsymbol{r})d\boldsymbol{r}
$$

を定義すると、最終的に固有エネルギー、そして係数$b_m$を求めるための方程式は以下の一般化固有値方程式で表され、方程式を解くことでBloch波数$\boldsymbol{k}$をパラメータとして展開に用いた関数の数だけ固有値・固有ベクトルが得られます。

$$
\begin{Vmatrix}
  M_{\boldsymbol{k}} - \varepsilon_{\boldsymbol{k}} S_{\boldsymbol{k}}
\end{Vmatrix} = 0
$$

・・・とはいえ、これではほぼ、ただの行列表示のシュレーディンガー方程式で具体的にどのような原子軌道を展開に用いればよいのかが分から無ければ計算もできません。

というわけで次節でいよいよ、Tight-bindingな描像に基づいた近似を導入していき、現実的に計算可能な方程式を導いていきます。前置きが長すぎ？

# 第二工程：Tight-binding近似

長かった前置きが終わり、ここから「電子が格子点の局所ポテンシャルにタイトにバインディングされている」描像のもとで以下のような近似を導入し、具体的な固有エネルギー、固有関数を求めていきます。



## Tight-bindingな近似その1：バンド$n$の展開に必要な原子軌道の見積もり^[この部分はアシュクラフト・マーミン上（I)を参考にした。]^[ものすごく余談ですが、最初にゼミで読んだ翻訳書がニジェ・オーランドでGrassmann代数をさらっとかじった後アルトランド・サイモンズで経路積分をかじってさらに例題を解いてくみたいな内容だったので、カナカナ連名著者の教科書は自分には敷居の高い高尚な内容が書かれているものだと思い込んでアシュクラフト・マーミン等も避けていたのですが、最近になって読んでみると非常に丁寧に知りたいことが書いてあり愕然としました。さっさと読めばよかった。。。]



まず初めに、重要な関係性として**各バンドに属する固有関数は少数の原子軌道関数で良く近似できる**ことを確かめます。^[私はこの部分が最近までイマイチ理解できておらず、Tight-bindingモデルの理解にかなり遠回りをしました]


さて、前節でWannier関数を原子軌道関数で（又はBloch関数をBloch和で）展開したのですが、展開に用いる原子軌道が多くなればなるほど精度が良くなるように思えます。それはそうなのですが、そんなことを言うと近似も近似でなくなるわけで、以下のようにある係数$b_m$、つまり展開に含まれる原子軌道のうち、$m$番目のエネルギー順位の原子軌道が含まれる割合についての関係式を得ることができ、結局そのエネルギーバンド程度の固有エネルギーに対応する（少数の）原子軌道関数のみで展開できることが分かります。

早速前に進みましょう。行列形式は一旦忘れて、初めに書いた変分法のもと、固有状態の試行関数が満たすべき方程式

$$
\left(
 -\frac{\hbar^2}{2m}\nabla^2
    +
   \sum_{\boldsymbol{R}}V(\boldsymbol{r} - \boldsymbol{R})
   \right)
   \varphi_{\boldsymbol{k}}(\boldsymbol{r}) = \varepsilon_{\boldsymbol{k}}\varphi_{\boldsymbol{k}}(\boldsymbol{r})  \\
   
\Rightarrow
\left(
 -\frac{\hbar^2}{2m}\nabla^2
    +
   \sum_{\boldsymbol{R}}V(\boldsymbol{r} - \boldsymbol{R})
   \right)
   \sum_{\boldsymbol{R}}e^{i\boldsymbol{k}\cdot\boldsymbol{R}} \sum_mb_m\phi_m(\boldsymbol{r}-\boldsymbol{R})
    = \varepsilon_{\boldsymbol{k}}\sum_{\boldsymbol{R}}e^{i\boldsymbol{k}\cdot\boldsymbol{R}} \sum_mb_m\phi_m(\boldsymbol{r}-\boldsymbol{R})
$$

に、今度は（Bloch和ではなく）原点を中心とした原子軌道関数の複素共役$\phi_{l}^*(\boldsymbol{r})$（$\boldsymbol{R} = \boldsymbol{0}$の場合に対応する）をかけて全空間について積分します。

まずは比較的シンプルな右辺について、右辺後の計算がわかりやすいように$\sum_{\boldsymbol{R}}$を$\boldsymbol{R} = \boldsymbol{0}$とそれ以外に分けて書くと、

$$
\begin{align*}
（右辺の積分）&=\varepsilon_{\boldsymbol{k}}
\left(
   \sum_m b_m \int \phi_{l}^*(\boldsymbol{r})\phi_{m}(\boldsymbol{r})d\boldsymbol{r}

   +
   \sum_{\boldsymbol{R}\neq\boldsymbol{0}}\sum_m b_m e^{i\boldsymbol{k}\cdot\boldsymbol{R}}
   
   \int \phi_{l}^*(\boldsymbol{r})\phi_{m}(\boldsymbol{r} - \boldsymbol{R})d\boldsymbol{r}

\right)\\

&=
\varepsilon_{\boldsymbol{k}}
\left(
    b_l
   +
   \sum_m 
   \left[
   \sum_{\boldsymbol{R}\neq\boldsymbol{0}}
   e^{i\boldsymbol{k}\cdot\boldsymbol{R}}
   
   \int \phi_{l}^*(\boldsymbol{r})\phi_{m}(\boldsymbol{r} - \boldsymbol{R})d\boldsymbol{r}
   \right]b_m 

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
d\boldsymbol{r}
\end{align*}
$$

特に$\varphi$が$\hat{H}$の固有状態のとき、固有値を$\varepsilon$とすると固有値は実数なので、$\int \left(\hat{H}\varphi(\boldsymbol{r})\right)^* \psi(\boldsymbol{r})d\boldsymbol{r} = \varepsilon\int \varphi^*(\boldsymbol{r}) \psi(\boldsymbol{r})d\boldsymbol{r}$、を利用して積分を計算していきます。

こちらも$\boldsymbol{R} = \boldsymbol{0}$の場合とそれ以外、さらにハミルトニアン部分も孤立原子のハミルトニアンとそれ以外で分けて

$$
\begin{align*}
(左辺の積分)&=
\sum_{\boldsymbol{R}}\sum_me^{i\boldsymbol{k}\cdot\boldsymbol{R}} b_m
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
      \right] b_m \\

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
      \right]  b_m 
      \\

&\>\>\>\>+
\sum_m
\left[
    \int \phi_l^*(\boldsymbol{r})
   \left(
   \sum_{\boldsymbol{R}'\neq \boldsymbol{0}}
      V(\boldsymbol{r}-\boldsymbol{R}')
      \right) \phi_m(\boldsymbol{r})d\boldsymbol{r}
      \right] b_m 
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
      \right] b_m\\

&=
\varepsilon_l^{\rm a}  b_l \\

&\>\>\>\>+
\varepsilon_l^{\rm a}
\sum_m 
\left[
   \sum_{\boldsymbol{R}\neq\boldsymbol{0}}
   e^{i\boldsymbol{k}\cdot\boldsymbol{R}} 
    \int \phi_l^*(\boldsymbol{r})
   \phi_m(\boldsymbol{r}-\boldsymbol{R})d\boldsymbol{r}
   \right] b_m\\



&\>\>\>\>+
\sum_m 
\left[
    \int \phi_l^*(\boldsymbol{r})
   \left(
   \sum_{\boldsymbol{R}'\neq \boldsymbol{0}}
      V(\boldsymbol{r}-\boldsymbol{R}')
      \right) \phi_m(\boldsymbol{r})d\boldsymbol{r}
      \right] b_m\\

&\>\>\>\>+
\sum_m 
\left[
   \sum_{\boldsymbol{R}\neq\boldsymbol{0}}
   e^{i\boldsymbol{k}\cdot\boldsymbol{R}} \int \phi_l^*(\boldsymbol{r})
   \left(
   \sum_{\boldsymbol{R}'\neq \boldsymbol{0}}
      V(\boldsymbol{r}-\boldsymbol{R}')
      \right)  \phi_m(\boldsymbol{r}-\boldsymbol{R})d\boldsymbol{r}
      \right]b_m  

\end{align*}
$$

となります。3式目の第1項、第2項の変形に、先ほど述べた原子軌道関数が孤立原子の固有状態であること、および固有値が実数であることを使いました。

さて、先ほどの等式から$b_l$単体の項とそれ以外をそれぞれ移項したりして整理すると、展開係数$b_l$に関する等式

$$
\begin{align*}
(\varepsilon_{\boldsymbol{k}} - \varepsilon_l^{\rm a})b_l 

&=
-(\varepsilon_{\boldsymbol{k}} - \varepsilon_l^{\rm a})
\sum_m 
\left[
   \sum_{\boldsymbol{R}\neq\boldsymbol{0}}
   e^{i\boldsymbol{k}\cdot\boldsymbol{R}} 
   \int \phi_l^*(\boldsymbol{r})
   \phi_m(\boldsymbol{r}-\boldsymbol{R})d\boldsymbol{r}
   \right]b_m \\

&\>\>\>\>+
\sum_m 
\left[
    \int \phi_l^*(\boldsymbol{r})
   \left(
   \sum_{\boldsymbol{R}'\neq \boldsymbol{0}}
      V(\boldsymbol{r}-\boldsymbol{R}')
      \right) \phi_m(\boldsymbol{r})d\boldsymbol{r}
      \right] b_m\\

&\>\>\>\>+
\sum_m 
\left[
   \sum_{\boldsymbol{R}\neq\boldsymbol{0}}
   e^{i\boldsymbol{k}\cdot\boldsymbol{R}} \int \phi_l^*(\boldsymbol{r})
   \left(
   \sum_{\boldsymbol{R}'\neq \boldsymbol{0}}
      V(\boldsymbol{r}-\boldsymbol{R}')
      \right)  \phi_m(\boldsymbol{r}-\boldsymbol{R})d\boldsymbol{r}
      \right]b_m 
\end{align*}
$$

を得ることができます。

等式の右辺についてひとつづつ見ていくと、まず第1項は異なる場所に中心を持つ原子軌道関数同士の積分が含まれています。ここで今、「格子点に強く束縛された状態」を想定していたことから、その展開に含まれる原子軌道関数も各格子点に局在していると考え、この積分の値は小さいとみなせます。
なお、ここで現れる積分

$$
\int \phi_l^*(\boldsymbol{r})
   \phi_m(\boldsymbol{r}-\boldsymbol{R})d\boldsymbol{r}
$$

は **Overlap Integral（重なり積分）** と呼ばれます。

また第2項についても、原点を中心とする原子軌道関数$\phi(\boldsymbol{r})$は格子間距離$\boldsymbol{R}$離れた場所では小さくなる一方、原点においては原点以外を中心として持つ孤立原子ポテンシャル$V(\boldsymbol{r}-\boldsymbol{R})$が小さくなると考えられるので、こちらの積分も小さい値を持つと考えられます。
^[アシュクラフト・マーミンでは「孤立原子ポテンシャルは波動関数と比べるとそこまで急速に小さくなるわけではないのでこの項が「小さい」とするのは幾分怪しいが、$\boldsymbol{k}$に依存しない項なので問題ないのと、ポテンシャルと波動関数を適切に選べばこの項も小さくできる」というようなことが脚注で書いてある。$\boldsymbol{k}$に依存しないとなぜ問題ないのか、どうやって小さくするのか等まだよくわかっていないのだがひとまず受け入れておくことにする]

最後に第3項についても、第1項、第2項と似たような理由で、中心が異なる（原点、$\boldsymbol{R}'$、$\boldsymbol{R}$）3つの局在した関数の積分となるので値は小さいと考えられます（これを「3中心積分などと呼んだりする」）。

以上のように、右辺はどの項も小さい値を持つと見積もることができ、左辺

$$
(\varepsilon_{\boldsymbol{k}} - \varepsilon_l^{\rm a})b_l
$$

は小さくならなければならないことがわかります。これを満たすには

$$
\varepsilon_{\boldsymbol{k}} \simeq \varepsilon_l^{\rm a}
$$

すなわち固有値$\varepsilon_{\boldsymbol{k}}$に対して、係数$b_l$に対応した原子軌道関数の固有エネルギー$\varepsilon_l^{\rm a}$が近い値を持つ場合、または

$$
b_l \simeq 0
$$

の場合であると言えます。

[固体（結晶）中の電子状態](https://zenn.dev/ponzumai/articles/tight-binding-model-electrons-in-solids)の章で見たように、バンド指標$n$に対応する固体内の電子の固有値は、
- バンド指標$n$ごとに離散的な、すなわち程度離れたエネルギーを持ち
- エネルギー$\varepsilon_{n,\boldsymbol{k}}$は連続、かつ限定的な領域「ブリルアンゾーン」を周期として持つ周期関数

なのでした。つまり上記固有値$\varepsilon_{\boldsymbol{k}}$は、どのバンド指標で表されるにせよ、BZ内のBloch波数$\boldsymbol{k}$の関数として最大値と最小値を持ち、一定の幅に収まっているわけです。

従って、「あるバンド$n$のBloch関数$\varphi_{n,\boldsymbol{k}}(\boldsymbol{r})$またはWannier関数$w_{n,\boldsymbol{R}}(\boldsymbol{r})$の展開に用いる原子軌道関数$\phi_l(\boldsymbol{r})$」は、そのBloch関数が属するエネルギーバンドの幅と近い原子順位$\varepsilon_m^{\rm a}$を持つ原子軌道関数$\phi_m$のみで良く近似される、ということを意味しています！^[ここの表現がだいぶまどろっこしいのは多分私が完全に理解していないからのようなきがするので、いつかスッキリ書き直したい]
（逆に、この近似の下ではある原子順位の重ね合わせで表現されるBloch関数のバンドは、概ねその原子準位に近い範囲で値を持つと言えます。）

もう少し具体的に書くと、縮退した、又は近い原子準位を持つ原子軌道関数を$m = m_1 , m_2,  m_3$として、そのグループをひとまず$q$でラベルしておき、原子準位を$\varepsilon^{{\rm a}_q}_{m_i}$と置くと、原子準位$\varepsilon^{{\rm a}_q}_{m_i}$に近い固有値に対応するBloch関数$\varphi^q_{\boldsymbol{k}}$は、

$$
\varphi_{n_q,\boldsymbol{k}}(\boldsymbol{r}) \simeq \sum_{m = m_1,m_2,m_3}b_m^{n_q}
\sum_{\boldsymbol{R}}e^{i\boldsymbol{k}\cdot\boldsymbol{R}} \phi_m(\boldsymbol{r}-\boldsymbol{R})
$$

と近似できることを意味します。これは大きな進歩です。

### 余談その3
~~私のように~~余計なことが気になってしまうような人は、前掲の関係式

$$
\begin{align*}
(\varepsilon_{\boldsymbol{k}} - \varepsilon_l^{\rm a})b_l

&=
-(\varepsilon_{\boldsymbol{k}} - \varepsilon_l^{\rm a})
\sum_m 
\left[
   \sum_{\boldsymbol{R}\neq\boldsymbol{0}}
   e^{i\boldsymbol{k}\cdot\boldsymbol{R}} 
   \int \phi_l^*(\boldsymbol{r})
   \phi_m(\boldsymbol{r}-\boldsymbol{R})d\boldsymbol{r}
   \right]b_m \\

&\>\>\>\>+
\sum_m 
\left[
    \int \phi_l^*(\boldsymbol{r})
   \left(
   \sum_{\boldsymbol{R}'\neq \boldsymbol{0}}
      V(\boldsymbol{r}-\boldsymbol{R}')
      \right) \phi_m(\boldsymbol{r})d\boldsymbol{r}
      \right] b_m\\

&\>\>\>\>+
\sum_m 
\left[
   \sum_{\boldsymbol{R}\neq\boldsymbol{0}}
   e^{i\boldsymbol{k}\cdot\boldsymbol{R}} \int \phi_l^*(\boldsymbol{r})
   \left(
   \sum_{\boldsymbol{R}'\neq \boldsymbol{0}}
      V(\boldsymbol{r}-\boldsymbol{R}')
      \right)  \phi_m(\boldsymbol{r}-\boldsymbol{R})d\boldsymbol{r}
      \right]b_m 
\end{align*}
$$

第1項について、「とはいえ、エネルギー順位が高い原子軌道関数は空間的にも広がっているだろうから、そういう場合を想定すれば第1項は大きくなり、「右辺が全部小さいから左辺も小さい、従ってエネルギーバンド近傍の原子順位の軌道だけで上手く近似できる」という結論にはならないのではないか？」等と~~いらんことを~~考えてしまうかもしれません。

そんな高エネルギー状態を展開に含めるのは物理的にはおかしいし、そもそも「局在した」関数で近似しよう、というのが出発点だとしたら上のようなことが気になるのもおかしいのかもしれませんが、一応安心するために、以下のような怪しい議論をしてみることにします。
**なお、以下のようなことを書いてある文献は私は寡聞にして知らず、まるで出鱈目なことを書いている可能性があるので注意してください**

（この後の章^[追記予定]で、さらに怪しい議論を展開して重なり積分が小さいことを納得していく予定です）


まず、式の簡略化のため、右辺第1項の「準位$l$の原子軌道関数$\phi_l$の重なり積分（と係数の積）の総和」を
$$
S_l\equiv\sum_m 
\left[
   \sum_{\boldsymbol{R}\neq\boldsymbol{0}}
   e^{i\boldsymbol{k}\cdot\boldsymbol{R}} 
   \int \phi_l^*(\boldsymbol{r})
   \phi_m(\boldsymbol{r}-\boldsymbol{R})d\boldsymbol{r}
   \right]b_m 
$$

第2項、第3項を合わせた「原点以外の周期ポテンシャルを間に挟んだ積分（と係数の積）の相和」を

$$
V_l\equiv

\sum_m 
\left[
   \sum_{\boldsymbol{R}}
   e^{i\boldsymbol{k}\cdot\boldsymbol{R}} \int \phi_l^*(\boldsymbol{r})
   \left(
   \sum_{\boldsymbol{R}'\neq \boldsymbol{0}}
      V(\boldsymbol{r}-\boldsymbol{R}')
      \right)  \phi_m(\boldsymbol{r}-\boldsymbol{R})d\boldsymbol{r}
      \right]b_m 
$$

と置きます。

そして、バンド$n$と近い原子準位を持ち、したがって局在した（重なり積分が小さい）準位を$l$、バンドと離れ、大きい原子準位を持ち、したがって広がりの大きい（重なり積分の大きい）準位を$L$と置いて、展開係数$b_l$と$b_L$の比を求めてみます。


$$
\begin{align*}
\frac{b_L}{b_l} &= 
\frac{-S_L + (\varepsilon_{n,k} - \varepsilon_L^{\rm a})^{-1}V_L}{-S_l + (\varepsilon_{k} - \varepsilon_l^{\rm a})^{-1}V_l}\\

\end{align*}
$$

ここで、準位$l$の固有関数は局在性が強いと仮定したので、$S_l\simeq 0$とすると、

$$
\begin{align*}
\frac{b_L}{b_l} &\simeq
\frac{-S_L + (\varepsilon_{k} - \varepsilon_L^{\rm a})^{-1}V_L}{(\varepsilon_{k} - \varepsilon_l^{\rm a})^{-1}V_l}\\

&=
(\varepsilon_{k} - \varepsilon_l^{\rm a})\frac{-S_L }{V_l}
+\frac{\varepsilon_{k} - \varepsilon_l^{\rm a}}{\varepsilon_{k} - \varepsilon_L^{\rm a}} V_L

\end{align*}
$$

となります。ここで第2項も$\varepsilon_{k} \simeq \varepsilon_l^{\rm a}$と$\varepsilon_{k} \ll \varepsilon_L^{\rm a}$の仮定から小さくなり、無視して良さそうなので落として、最終的に「広がった」原子準位と「局在した」原子準位の比

$$
\begin{align*}
\frac{b^n_L}{b^n_l} &\simeq

(\varepsilon_{k} - \varepsilon_l^{\rm a})\frac{-S_L }{V_l}

\end{align*}
$$

を得ます。この式を眺めてみると、例え重なり$S_L$が大きくても、周期ポテンシャルを含む積分部分


$$
V_l=

\sum_m 
\left[
   \sum_{\boldsymbol{R}}
   e^{i\boldsymbol{k}\cdot\boldsymbol{R}} \int \phi_l^*(\boldsymbol{r})
   \left(
   \sum_{\boldsymbol{R}'\neq \boldsymbol{0}}
      V(\boldsymbol{r}-\boldsymbol{R}')
      \right)  \phi_m(\boldsymbol{r}-\boldsymbol{R})d\boldsymbol{r}
      \right]b_m 
$$

がそれなりに大きければ、

$$
\begin{align*}
\frac{b_L}{b_l} &\simeq

(\varepsilon_{k} - \varepsilon_l^{\rm a})\frac{-S_L }{V_l}
\simeq 0

\end{align*}
$$

とでき、先ほどの結論に問題はなさそうです。
これは局所ポテンシャル（孤立原子のポテンシャル）が大きく、空間的に離れた原子軌道関数の減衰を補う位の積分値を与えてくれるような場合、という描像に対応してそうです。
「固体内の電子が格子点に強く束縛されている＝孤立原子のポテンシャルが強い」という、当初の「Tight-binding」な仮定とも一致しており、良い感じです。
（先ほどはこの部分を「小さい」と考えて同様の結論を導いたので矛盾しているように思えますが、比べている対象が違っていて、先ほどはエネルギー差$\varepsilon_{\boldsymbol{k}} - \varepsilon_l^{\rm a}$と比べて「小さい」としていて、今回は（広がった原子軌道の）重なり積分の相和$S_L$と比べていますが、展開係数$b_l^n$についての規格化条件を考えると、$S_L$は最大でも$1$ですので、それに比べて（それなりに大きい）という条件はそこまで矛盾しないかと思います。

さて、逆に周期ポテンシャルが弱い場合、極端な場合では重なり積分の和$S_l$とほぼ一致する場合は、$S_L/V_l$の分子が小さくなり、$(\varepsilon_{k} - \varepsilon_l^{\rm a})$と打ち消しあっちゃって$b_L/b_l$が大きくなり得て、近似が上手くいかなさそうだと考えられます。

というわけで、それなりに強い局所ポテンシャルを考えている限りにおいては、安心してあるバンド$n$のBloch関数を少数のBloch和で（Wannier関数を少数の原子軌道関数で）展開できそうです。（まあそもそも、そこまでややこしいこと考えずに、展開してみて実験と合えばそれが正義なのかもしれませんが）

### 固有値方程式のブロック対角化

以上のようにWannier関数が少数の原子軌道で、つまりBloch関数が少数のBloch和で展開されると仮定した場合、Bloch関数$\varphi_{\boldsymbol{k}}(\boldsymbol{r})$の展開はエネルギー的に近い原子準位の集まりを用いてグループ分けできることがわかりました。
イメージ的には水素原子の固有関数のラベルを用いて、縮退した原子軌道を1グループとすると、

- $1s$グループ：$\varphi_{1s,\boldsymbol{k}}\simeq\sum_{\boldsymbol{R}}e^{i\boldsymbol{k}\cdot\boldsymbol{R}} \phi_{1s}(\boldsymbol{r}-\boldsymbol{R})$
- $2s$グループ：$\varphi_{2s,\boldsymbol{k}}\simeq\sum_{\boldsymbol{R}}e^{i\boldsymbol{k}\cdot\boldsymbol{R}} \phi_{2s}(\boldsymbol{r}-\boldsymbol{R})$
- $2p$グループ：$\varphi_{2p,\boldsymbol{k}}\simeq\sum_{m = 1,2,3}b_m^{2p}\sum_{\boldsymbol{R}}e^{i\boldsymbol{k}\cdot\boldsymbol{R}} \phi_{2p}(\boldsymbol{r}-\boldsymbol{R})$
- $3s$グループ：$\cdots$


となります。（なお、この段階では異なる軌道間の重なり積分がまだ値を持つので、異なるグループに属する、つまり異なるバンドに属する（近似）ブロッホ関数は直交しません。後に重なり積分をゼロと置く近似をすることで異なるバンド間のBloch関数の正規直交性が満たされることになります）

これにより、先ほど求めた固有値方程式もサイズの小さい行列のみを扱えばよいことになります。

すなわち、例えば$m = m_1, m_2, m_3$の原子準位に対応する原子軌道（このグループを先ほどと同様に$q$と名付けます）で展開した場合、つまり$\varepsilon_{n_q,\boldsymbol{k}} \simeq \varepsilon_{m_1}^{\rm a}, \varepsilon_{m_2}^{\rm a}, \varepsilon_{m_3}^{\rm a}$であるとき、Bloch関数の展開は

$$
\varphi_{n_q,\boldsymbol{k}}\simeq
\sum_{m = m_1, m_2, m_3}b_m^q
\left[\sum_{\boldsymbol{R}}e^{i\boldsymbol{k}\cdot\boldsymbol{R}} \phi_m(\boldsymbol{r}-\boldsymbol{R})
\right]
$$

となるため、工程1で求めた固有値方程式は対応する原子準位$\phi_{m_i}(\boldsymbol{r}), i = 1,2,3$を用いて、次の$3\times 3$

$$
(M_{\boldsymbol{k}}^q)_{ij} \equiv

$$

行列

$$
\left(M_{\boldsymbol{k}}^q\right)_{ij} \equiv \int \Phi_{m_i,\boldsymbol{k}}^{q*}(\boldsymbol{r})\hat{H}\Phi^{q}_{m_j,\boldsymbol{k}}(\boldsymbol{r})d\boldsymbol{r},\\
\left(S_{\boldsymbol{k}}^q\right)_{ij} \equiv \int \Phi_{n,\boldsymbol{k}}^*(\boldsymbol{r})\Phi_{m,\boldsymbol{k}}(\boldsymbol{r})d\boldsymbol{r}
$$

を定義すると、最終的に固有エネルギー、そして係数$b_m$を求めるための方程式は以下の固有値方程式で表され、方程式を解くことでBloch波数$\boldsymbol{k}$をパラメータとして展開に用いた関数の数だけ固有値・固有ベクトルが得られます。

$$
\begin{Vmatrix}
  M_{\boldsymbol{k}} - \varepsilon_{\boldsymbol{k}} S_{\boldsymbol{k}}
\end{Vmatrix} = 0
$$


$$
\begin{vmatrix}
(M_{\boldsymbol{k}})_{11} -\varepsilon_{n_q,\boldsymbol{k}}(S_{\boldsymbol{k}})_{11} & (M_{\boldsymbol{k}})_{12} -\varepsilon_{n_q,\boldsymbol{k}}(S_{\boldsymbol{k}})_{12} & (M_{\boldsymbol{k}})_{13} -\varepsilon_{n_q,\boldsymbol{k}}(S_{\boldsymbol{k}})_{13} \\

(M_{\boldsymbol{k}})_{21} -\varepsilon_{n_q,\boldsymbol{k}}(S_{\boldsymbol{k}})_{21} & (M_{\boldsymbol{k}})_{22} -\varepsilon_{n_q,\boldsymbol{k}}(S_{\boldsymbol{k}})_{22} & (M_{\boldsymbol{k}})_{23} -\varepsilon_{n_q,\boldsymbol{k}}(S_{\boldsymbol{k}})_{23} \\

(M_{\boldsymbol{k}})_{31} -\varepsilon_{n_q,\boldsymbol{k}}(S_{\boldsymbol{k}})_{31} & (M_{\boldsymbol{k}})_{32} -\varepsilon_{n_q,\boldsymbol{k}}(S_{\boldsymbol{k}})_{32} & (M_{\boldsymbol{k}})_{33} -\varepsilon_{n_q,\boldsymbol{k}}(S_{\boldsymbol{k}})_{33} \\
\end{vmatrix}
=0
$$

を解けばよいことになります。凄い省エネです。

上記のように、求めたいバンドに応じてエネルギー的に近そうな原子準位を選べば、その数程度のサイズの方程式（永年方程式）を解けばよいことになります。
とはいえ近似なので、完全に縮退した原子準位のみではなく、近いエネルギー順位も取り入れたりする必要はあるみたいですが、それにしても高々$4$個（$s+p$）とか、$6$個（$s+d$）位なので、十分計算できる範囲ですね。

これだけでもかなりの成果ですが、具体的な解を得るためにここからもう少し近似を入れていきます。



## Tight-bindingな近似その2：原子軌道が強く束縛されている感じの近似

さて、ここでより「強く束縛された状態」の描像に基づいた近似、つまり原子準位の広がりが小さいという近似を行っていきます。これにより簡単な場合であれば具体的なエネルギーバンド・固有関数が得られます。

### 重なり積分の近似



まず、先ほども考えたように、原子軌道$varphi_m(\boldsymbol{r}) $に対して、異なる格子点を中心に持つ原子軌道同士の「重なり」が小さいとして、思い切って次の積分を$0$と考えます。つまり


$$
\int\varphi_m^*(\boldsymbol{r}-\boldsymbol{R}) \varphi_l(\boldsymbol{r-\boldsymbol{R}'})d\boldsymbol{r} \simeq \delta_{m,n}\delta_{\boldsymbol{R},\boldsymbol{R}'}
$$

とします。この積分（Overlap Integral、重なり積分）が小さいというのは格子間距離に比べてWannier関数を近似している原子軌道の広がりが十分小さい、つまり「強く束縛されている」ようなイメージに対応しています。

これによりLCAO近似で表したWannier関数も正規直交性を持ち、Bloch関数も正規直交性を満たすようになります。

そして、先ほど得た永年方程式において、行列$S_{\boldsymbol{k}}$が単位行列になる：

$$
\begin{align*}
\left(S_{\boldsymbol{k}}\right)_{nm} &= \int \Phi_{n,\boldsymbol{k}}^*(\boldsymbol{r})\Phi_{m,\boldsymbol{k}}(\boldsymbol{r})d\boldsymbol{r}
\\

&=

\sum_{\boldsymbol{R}'}
\sum_{\boldsymbol{R}}e^{i\boldsymbol{k}\cdot(\boldsymbol{R}'-\boldsymbol{R})}\int
 \phi_n^*(\boldsymbol{r}-\boldsymbol{R}')
 \phi_m(\boldsymbol{r}-\boldsymbol{R})d\boldsymbol{r}\\

 &=
 \sum_{\boldsymbol{R}'}
\sum_{\boldsymbol{R}}e^{i\boldsymbol{k}\cdot(\boldsymbol{R}'-\boldsymbol{R})}
\delta_{n,m}\delta_{\boldsymbol{R}',\boldsymbol{R}}\\

&=
N\delta_{n,m}

\end{align*}
$$

ことで、永年方程式は

$$
\begin{Vmatrix}
(M_{\boldsymbol{k}})_{nm} - N\varepsilon_{n,\boldsymbol{k}}\delta_{nm}
\end{Vmatrix}
=0
$$

となり、行列$M_{\boldsymbol{k}}/N$の対角化に帰着します。

### 飛び移り積分

さらに

$$
\begin{align*}
\left(M_{\boldsymbol{k}}\right)_{nm} 
&=
 \int \Phi_{n,\boldsymbol{k}}^*(\boldsymbol{r})\hat{H}\Phi_{m,\boldsymbol{k}}(\boldsymbol{r})d\boldsymbol{r},\\

&=
\sum_{\boldsymbol{R}_i}
\sum_{\boldsymbol{R}_j}e^{i\boldsymbol{k}\cdot(\boldsymbol{R}_i-\boldsymbol{R}_j)}\int
 \phi_n^*(\boldsymbol{r}-\boldsymbol{R}_i)
 \left(
 -\frac{\hbar^2}{2m}\nabla^2
    +
   \sum_{\boldsymbol{R}_l}V(\boldsymbol{r} - \boldsymbol{R}_l)
   \right)
 \phi_m(\boldsymbol{r}-\boldsymbol{R}_j)d\boldsymbol{r}\\
\end{align*}
$$

とできますが、ここで並進の対称性を利用して$N$個の$\boldsymbol{R}_i$に対して積分は同じ値を取るので、$\boldsymbol{R}_i = \boldsymbol{0}$の場合の$N$倍を考えることで

$$
\begin{align*}
&=
N
\sum_{\boldsymbol{R}_j}e^{i\boldsymbol{k}\cdot(\boldsymbol{R}_j)}\int
 \phi_n^*(\boldsymbol{r})
 \left(
 -\frac{\hbar^2}{2m}\nabla^2
    +
   \sum_{\boldsymbol{R}_l}V(\boldsymbol{r} - \boldsymbol{R}_l)
   \right)
 \phi_m(\boldsymbol{r}-\boldsymbol{R}_j)d\boldsymbol{r}\\

\end{align*}
$$

となります。さらに原子軌道関数が孤立原子のハミルトニアンの固有関数であることを利用し、また$\boldsymbol{R}_l = \boldsymbol{R}_j$とそれ以外の場合で分けることで

$$
\begin{align*}

 &=
N
\sum_{\boldsymbol{R}_j}e^{i\boldsymbol{k}\cdot\boldsymbol{R}_j}
\varepsilon_{n}^{\rm a}
\int
 \phi_n^*(\boldsymbol{r})
 
 \phi_m(\boldsymbol{r}-\boldsymbol{R}_j)d\boldsymbol{r}\\

&\>\>\>\>+
N
\sum_{\boldsymbol{R}_j}e^{i\boldsymbol{k}\cdot\boldsymbol{R}_j}

\int
 \phi_n^*(\boldsymbol{r})
 \left(
   \sum_{\boldsymbol{R}_l\neq \boldsymbol{0}}V(\boldsymbol{r} - \boldsymbol{R}_l)
   \right)
 \phi_m(\boldsymbol{r}-\boldsymbol{R}_j)d\boldsymbol{r}\\

 &=
N
\varepsilon_{n}^{\rm a}
\delta_{n,m}\\

&\>\>\>\>+
N

\int
 \phi_n^*(\boldsymbol{r})
 \left(
   \sum_{\boldsymbol{R}_l\neq \boldsymbol{0}}V(\boldsymbol{r} - \boldsymbol{R}_l)
   \right)
 \phi_m(\boldsymbol{r})d\boldsymbol{r}\\

&\>\>\>\>+
N
\sum_{\boldsymbol{R}_j}
e^{i\boldsymbol{k}\cdot\boldsymbol{R}_j}
\int
 \phi_n^*(\boldsymbol{r})
 
   V(\boldsymbol{r} - \boldsymbol{R}_j)
 
 \phi_m(\boldsymbol{r}-\boldsymbol{R}_j)d\boldsymbol{r}\\

&\>\>\>\>+
N
\sum_{\boldsymbol{R}_j}
e^{i\boldsymbol{k}\cdot\boldsymbol{R}_j}
\int
 \phi_n^*(\boldsymbol{r})
 \left(
   \sum_{\boldsymbol{R}_l\neq \boldsymbol{0}, \boldsymbol{R}_j}V(\boldsymbol{r} - \boldsymbol{R}_l)
   \right)
 \phi_m(\boldsymbol{r}-\boldsymbol{R}_j)d\boldsymbol{r}\\
\end{align*}
$$

となります。

ここで、第2項の積分を

$$
\begin{align*}
\int
 \phi_n^*(\boldsymbol{r})
 \left(
   \sum_{\boldsymbol{R}_l\neq \boldsymbol{0}}V(\boldsymbol{r} - \boldsymbol{R}_l)
   \right)
 \phi_m(\boldsymbol{r})d\boldsymbol{r}
 
 &\equiv
\Delta\varepsilon_{nm}
 \end{align*}
$$

と定義します。この積分は（グロッソ・パラビチニによると）「結晶場積分」と呼ばれているそうなのですが、なんとググっても1件もヒットしないので、あまり一般的な用語ではないのかもしれません。（"Crystal Fiels Integral(s)"では1000件ほど（のみ）ヒット。）

なお、原点以外に中心を持つ局所ポテンシャルと、原点を中心とする原子軌道の内積の積を積分していますが、「中心からそれなりに離れた場所では原子ポテンシャルは大体一定になっているだろう」という近似のもと、

$$
\begin{align*}
\int
 \phi_n^*(\boldsymbol{r})
 \left(
   \sum_{\boldsymbol{R}_l\neq \boldsymbol{0}}V(\boldsymbol{r} - \boldsymbol{R}_l)
   \right)
 \phi_m(\boldsymbol{r})d\boldsymbol{r}
 
 &\simeq
\Delta\varepsilon\delta_{nm}
 \end{align*}
$$

と定数の対角行列と近似されることもあるようです。このようにすればこの項はエネルギーを定数シフトさせるだけの効果となり、後の計算では無視していけるようになります。（本章では一応残しておきます。）

第3項は、原点と異なる格子点と、それと等しい局所ポテンシャルの積で、関数の局在性から近い位置の間の積分のみ値を持つと考え、特に原点から**再隣接格子点**へのベクトルを$\boldsymbol{R}_I$として

$$
\sum_{\boldsymbol{R}_j}
e^{i\boldsymbol{k}\cdot\boldsymbol{R}_j}
\int
 \phi_n^*(\boldsymbol{r})
 
   V(\boldsymbol{r} - \boldsymbol{R}_j)
 
 \phi_m(\boldsymbol{r}-\boldsymbol{R}_j)d\boldsymbol{r}

 \simeq
\sum_{\boldsymbol{R}_I}
e^{i\boldsymbol{k}\cdot\boldsymbol{R}_I}
\int
 \phi_n^*(\boldsymbol{r})
 
   V(\boldsymbol{r} - \boldsymbol{R}_I)
 
 \phi_m(\boldsymbol{r}-\boldsymbol{R}_I)d\boldsymbol{r}

$$

と近似されます。ここで$\sum_{\boldsymbol{R}_I}$は、再隣接格子点のみ和を取ることを表します。

積分

$$
\int
 \phi_n^*(\boldsymbol{r})
 
   V(\boldsymbol{r} - \boldsymbol{R}_I)
 
 \phi_m(\boldsymbol{r}-\boldsymbol{R}_I)d\boldsymbol{r}
 \equiv
 t_I^{n,m}
$$

は、上記のように$t$で表示されることが多く、 **"Transfer Integral"、"Hopping Integral"、「飛び移り積分」** 等と呼ばれます。（多分"transfer"の"t"だと思います）
多くの場合再隣接格子、時には次近接辺りまでの値が近似で採用されます。


また第4項は、全て中心が異なる関数の積分であり、同様にそれぞれの関数が局在していることから無視することができると考えます。（この形の積分を「3中心積分」等と呼ぶようです）

以上をまとめると、対角化すべき行列$M_{\boldsymbol{k}}/N$は

$$
\frac{M_{\boldsymbol{k}}}{N} =

\begin{bmatrix}
\varepsilon_{m_1}^{\rm a} - \Delta\varepsilon_{m_1m_1} & \sum_{\boldsymbol{R}_I}e^{i\boldsymbol{k}\cdot\boldsymbol{R}_I} t_I^{m_1m_2} -\Delta\varepsilon_{m_1m_2} & \sum_{\boldsymbol{R}_I}e^{i\boldsymbol{k}\cdot\boldsymbol{R}_I} t_I^{m_1m_3}- \Delta\varepsilon_{m_1m_3}\\

\sum_{\boldsymbol{R}_I}e^{i\boldsymbol{k}\cdot\boldsymbol{R}_I} t_I^{m_2m_1} - \Delta\varepsilon_{m_2m_1} &
  \varepsilon_{m_2}^{\rm a} -\Delta\varepsilon_{m_2m_2} &
   \sum_{\boldsymbol{R}_I}e^{i\boldsymbol{k}\cdot\boldsymbol{R}_I} t_I^{m_2m_3} - \Delta\varepsilon_{m_2m_3}\\

\sum_{\boldsymbol{R}_I}e^{i\boldsymbol{k}\cdot\boldsymbol{R}_I} t_I^{m_3m_1} - \Delta\varepsilon_{m_3m_1} &
 \sum_{\boldsymbol{R}_I}e^{i\boldsymbol{k}\cdot\boldsymbol{R}_I} t_I^{m_3m_2} - \Delta\varepsilon_{m_3m_2} & 
 \varepsilon_{m_3}^{\rm a} - \Delta\varepsilon_{m_3m_3}
\end{bmatrix}
$$

となります。局所ポテンシャル$V(\boldsymbol{r})$や原子軌道関数$\phi_m(\boldsymbol{r})$は既に得られているものとしているので、行列の各要素は（少なくとも）数値的に得られて、後はサイズの小さい行列を対角化することで、展開に用いたエネルギー的に近い原子軌道の数（今回は3個）のエネルギーバンドと、3つの原子軌道の重ね合わせで表現される固有関数を得ることができます！

### 原子軌道や分子軌道の名前を用いたバンドの呼び方

この時、例えば展開に$s$軌道を選んだ場合、そこから得られるバンドを「$s$バンド」や「$s$的なバンド」、エネルギーが縮退した3つの$p$軌道で展開した場合に得られるバンドを「$p$バンド」「$p$的なバンド」などと呼びます。
また、（私はあまり詳しくないですが）展開の種類によっては分子軌道の用語「$\sigma$バンド」や「$\pi$バンド」なんかもよく見かけます。ここで「$\sigma$バンド」は$s, p_x, p_y$の展開、「$\pi$バンド」は残りの$p_z$軌道による展開に対応しているようです。


# Tight-binding Model

さて、上記のようにして固有値を求める方程式を得ることができました。この時、結晶場積分の値が全て等しいとして$\Delta\varepsilon\delta_{nm}$とおき、かつ展開に用いた軌道関数はエネルギーが近いものでしたので、思い切って原子準位$\varepsilon_{m_i}^{\rm a}$も全て等しいと近似し、その値を$\varepsilon^{\rm a}$と置けば、単位行列$I$を用いて

$$
\frac{M_{\boldsymbol{k}}}{N} = (\varepsilon^{\rm a} - \Delta\varepsilon)I + 
\begin{bmatrix}
0 & \sum_{\boldsymbol{R}_I}e^{i\boldsymbol{k}\cdot\boldsymbol{R}_I} t_I^{m_1m_2}  &
 \sum_{\boldsymbol{R}_I}e^{i\boldsymbol{k}\cdot\boldsymbol{R}_I} t_I^{m_1m_3}\\

\sum_{\boldsymbol{R}_I}e^{i\boldsymbol{k}\cdot\boldsymbol{R}_I} t_I^{m_2m_1} &
  0 &
   \sum_{\boldsymbol{R}_I}e^{i\boldsymbol{k}\cdot\boldsymbol{R}_I} t_I^{m_2m_3} \\

\sum_{\boldsymbol{R}_I}e^{i\boldsymbol{k}\cdot\boldsymbol{R}_I} t_I^{m_3m_1}  &
 \sum_{\boldsymbol{R}_I}e^{i\boldsymbol{k}\cdot\boldsymbol{R}_I} t_I^{m_3m_2}  & 
 0
\end{bmatrix}
$$

となり、エネルギーバンドの形状は飛び移り積分$t_I^{nm}$の値（と結晶構造）のみで決定されることになります。つまり飛び移り積分をパラメータとして見ることができるようになります。

そして、ここからは統計力学の範疇になるのですが、エネルギーバンドの形状が分かれば、様々な固体の物理量を計算することができるようになります。従ってエネルギーバンドは、固体の性質を端的に表現する関数であり、固体を表す「モデル」だともいえるわけで、上記のように考えるとそのモデルがパラメータとして$t$を用いて表現された形になってると見れます。

というわけで上記のような$t$をパラメータとして表したエネルギーバンドを指して、Tight-binding Modelということもあります。


## 簡単な例でのエネルギーバンド＝Tight-bindingモデル

最後に、最も簡単な例として、1次元格子において、固有関数を$s$軌道を用いて展開した場合について具体的に考えてみます。$s$軌道は縮退がないので、もっとも粗い近似では1つの原子軌道（Bloch和）のみで、Wannier関数（Bloch関数）を展開できます。

なお、1次元系をどのように表現すればよいかちょっとややこしいのですが、最も手っ取り早い方法は基本格子ベクトル$\boldsymbol{a}_i$を、

$$
|\boldsymbol{a}_2| = |\boldsymbol{a}_3| = A \gg |\boldsymbol{a}_1|
$$

と置いてしまうことかと思います。これは実際の1次元系について、それ以外の方向との格子間距離が離れている状況を表しています。（2方向に結晶サイズ1、つまり$N_2 = N_3 =1$の周期的境界条件を課しても上手くいきそうな気もするのですが、なんかややこしいのでやめておきます。）こうすると後述するように、考慮しない方向の飛び移り積分が$0$となり、1次元的なバンドが得られます。

さて、格子間隔が開いたとはいえ、固有関数（Bloch関数）は先述のように一つの$s$軌道関数を用いて、

$$
\varphi_{s,\boldsymbol{k}}(\boldsymbol{r}) \simeq \sum_{\boldsymbol{R}}e^{i\boldsymbol{k}\cdot\boldsymbol{R}}\phi_s(\boldsymbol{r}-\boldsymbol{R})
$$

と表せます。ここでややこしいですが、格子は1次元でも、その中の（固有）関数は3次元的に広がっており、展開に用いる原子軌道も3次元関数です。

1つの軌道で展開した場合は行列$M_{\boldsymbol{k}}$は行列ではなくなり、固有値方程式は


$$
\varepsilon_{s}^{\rm a} - \Delta\varepsilon_{s} + \sum_{\boldsymbol{R}_I}e^{i\boldsymbol{k}\cdot\boldsymbol{R}_I} t_I
 - \varepsilon_{s,\boldsymbol{k}}
=0
\\

\Rightarrow
\varepsilon_{s,\boldsymbol{k}}
=
\varepsilon_{s}^{\rm a} - \Delta\varepsilon_{s} + \sum_{\boldsymbol{R}_I}e^{i\boldsymbol{k}\cdot\boldsymbol{R}_I} t_I
$$

と、すでに対角化された等式になります。
ここで$\varepsilon_{s}^{\rm a}$は$s$軌道の原子準位、$\Delta\varepsilon_{s}$は$s$軌道の結晶場積分、そして$\varepsilon_{s,\boldsymbol{k}}$が、電子の固有状態が持つ固有エネルギーの式（エネルギーバンド）を表します。


ここでここまではまだ3次元っぽいですが、ここで飛び移り積分$t_I$について考えると、1次元系っぽくなります。
すなわち、$\boldsymbol{a}_2, \boldsymbol{a}_3$方向の再隣接格子へのベクトルを$\boldsymbol{A}$とし、その時の飛び移り積分を$t_A$と置くと、

$$
t_A = \int \phi_s^*(\boldsymbol{r})V(\boldsymbol{r} \pm \boldsymbol{A})\phi_s(\boldsymbol{r}\pm \boldsymbol{A})d\boldsymbol{r}
$$

となりますが、ここで$A$は大きい数と置いたことから、$\boldsymbol{a}_2, \boldsymbol{a}_3$方向の飛び移り積分は$t_A\simeq0$となります。この結果エネルギーバンドは

$$
\varepsilon_{s,\boldsymbol{k}}
\simeq
\varepsilon_{s}^{\rm a} - \Delta\varepsilon_{s} + \sum_{R_1 = \pm a}e^{ik_1 R_1} t_I
$$

と、波数空間における1次元関数（$k_2, k_3$方向には定数関数）のエネルギーバンドが得られます。

さらに、$V(\boldsymbol{r})$が球対称ポテンシャルを仮定しており、また$s$軌道も角度部分は定数、つまり球対称であるため、$R_I = \pm a$いずれも同じ積分の値を持ちます。
最後に、$s$軌道は実数関数かつ正、つまり$\phi_s(\boldsymbol{r})^* = \phi_s(\boldsymbol{r}) > 0$であることと、$V$は引力ポテンシャルであることから$V(\boldsymbol{r})<0$より、積分は負の値になります。

そこで

$$
-\int \phi_s^*(\boldsymbol{r})V(\boldsymbol{r} + \boldsymbol{a}_1)\phi_s(\boldsymbol{r} + \boldsymbol{a}_1)
=

-\int \phi_s^*(\boldsymbol{r})V(\boldsymbol{r} - \boldsymbol{a}_1)\phi_s(\boldsymbol{r} - \boldsymbol{a}_1)
\equiv t >0
$$

と置いて、1次元関数の固有エネルギー（エネルギーバンド）は

$$
\begin{align*}
\varepsilon_{s}(k)
&=
\varepsilon_{s}^{\rm a} -\Delta\varepsilon_{s}-t \sum_{R_I = \pm a}e^{ikR} \\

&=
\varepsilon_{s}^{\rm a} - \Delta\varepsilon_{s} -2t \cos(ka)
\end{align*}
$$

と、$s$軌道の原子順位$\varepsilon_s^{\rm a}$（から結晶場積分を引いたもの）を中心に、飛び移り積分$t$程度に広がったバンド構造をとることがわかります。



# おわりに

随分長くなってしまいましたが、本章では

- Bloch関数を原子軌道関数（Bloch和）で展開し、その際の固有値方程式を求め
- あるバンドに属するBloch関数はお互いに近い原子準位を持つ少数の原子軌道関数のみを用いて近似的に展開できることを示し
- 飛び移り積分を考えることでエネルギーバンドを


後は、本章で分かったこと（$+ \alpha$）をもとにしながら、第二量子化表示のTight-bindingモデルを導けば、ようやく本稿のゴールに到達です。

ただ、その前に少しだけ寄り道をして、次回・次々回ではTight-binding近似で定義された「飛び移り積分」「重なり積分」の物理的な意味を考えてみたいと思います。