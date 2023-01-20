---
title: "第一量子化のTight-bindingモデル（前編）"
emoji: "🧵"
type: "tech" # tech: 技術記事 / idea: アイデア
topics: ["quantum","quantumcomputing","quantumcomputer","物理","物理学"]
published: true
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
- 正規直交性：$\int_V w^*(\boldsymbol{r}-\boldsymbol{R}')w(\boldsymbol{r} - \boldsymbol{R})dr =\delta_{\boldsymbol{R},\boldsymbol{R}'}$
- 局在性：ラベル＝格子点の座標$\boldsymbol{R}$を中心として局在した関数

を満たすのでした。なお、ここではバンド指標$n$は指定せずに書いています。

## Tight-bindingモデルの構築について

これまでの振り返りを踏まえて、いくつかの仮定や近似を入れながらTight-bindingモデルを構築していきます。

その道のりは、大きく2つの工程に分かれます。

第一に、固体内のポテンシャルを孤立原子ポテンシャルで近似し、かつ固有関数であるBloch関数を原子軌道関数（孤立原子ポテンシャルの固有関数）の線形結合で表し固体内の電子の1体固有関数を求める方程式を導出します。

第二に、上記方程式に対して、「電子が孤立原子に強く束縛されている」と仮定した近似を行い、少数の原子軌道による展開（**LCAO近似**）のもと、具体的な固有エネルギー（エネルギーバンド）を導きます。

書いていくとかなり長くなってしまったので、本章では上記の第二工程前半LCAO近似までを説明します。

**（余談）** この記事を書いていて、色々参考にしていると、同じような操作でも「Tight-binding（又は強束縛、強結合等々）近似」と書いていたり、「Tight-binding（又は以下略）モデル」と書いていたり、どういう部分が「近似」と呼ばれてどこまで行くと「モデル」になるのかよくわからなくなってきました。
その辺の厳密な使い分けがあるのかどうか今のところよくわかっていないのですが、本章では、
- 「電子が強く束縛されている」という描像に基づいた近似を「Tight-binding近似」
- そのような近似によって得られた最終的な固有エネルギー（エネルギーバンド）を「Tight-bindingモデル」

と呼び分けることにします。全然違ってたら、すみません。


# 第一工程：周期ポテンシャルの近似＋原子軌道による展開$\rightarrow$固有値方程式の導出
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


## 原子軌道関数の線形結合を用いたWannier関数の展開（LCAO法） 

上記のようにして求めた固有関数形$\{\phi_m\}$は（多分）^[孤立原子シュレーディンガー方程式の固有関数系が完全系となることについてちょっと自信がないのですが、例えば動径関数を何らかの完全系$\{\varphi_n(r)\}$の線形結合で、角度部分は球面調和関数で展開したとして、そのように展開された固有関数系を全て（無限個）集めたとしたら、それは結局行列表示した孤立原子ハミルトニアンを対角化する行列による、完全系$\{\psi_l(\boldsymbol{r})\}( = \{\varphi_n(r) Y_l^m(\theta,\phi)\})$のユニタリ変換とみなせて、それは結局完全系を満たす、というようなことなのかと思っています。なんか出鱈目なことを書いていそうで怖いですが]
完全系を成すので、**固体中の電子の波動関数を原子軌道関数で展開**することができます。

特に、前章で求めたように、固体内の電子は「格子点に局在した関数」Wannier関数の線形結合で表現できるのでした。このWannier関数を、孤立原子（又は分子等）の波動関数$\phi_m(\boldsymbol{r})$と展開係数$b_m$を用いて、

$$
w_{\boldsymbol{R}}(\boldsymbol{r}) = w(\boldsymbol{r}-\boldsymbol{R}) = \sum_mb_m\phi_m(\boldsymbol{r}-\boldsymbol{R}) 
$$

と展開します。$b_m$は規格化条件$\sum_m|b_m|^2 = 1$を満たすものとします。
ここで、前章ではバンド指標を含めたBloch関数とそのWannier関数での展開を考えましたが、ここではバンドによらない形で定義しておきます。（後に固有値問題を解くことで各バンドに対応したBloch関数を求めます。

また、原子軌道関数であることがわかりやすいように、以降、原子軌道関数は先ほどと同じく$\phi$を用いて表します。

このような考え方は、「原子軌道の線形結合」の英語バージョン"Linear Combination of Atomic Orbitals"の頭文字をとってLCAO法、またはLCAO近似（展開する原子軌道に制限をつけなければ「近似」ではないですが、多くの場合物理的な考察に応じて少数の原子軌道で展開することになります。（そうしないとそもそも計算できないし）そのように少数の原子軌道で展開するときに「LCAO近似」というのだと思います）と呼ばれています。

この結果、結晶内の電子の固有状態＝Bloch関数は

$$
\begin{align*}
\varphi_{\boldsymbol{k}} &= \sum_{\boldsymbol{R}}e^{i\boldsymbol{k}\cdot\boldsymbol{R}}w(\boldsymbol{r}-\boldsymbol{R})\\

&=
 \sum_{\boldsymbol{R}}e^{i\boldsymbol{k}\cdot\boldsymbol{R}} \sum_mb_m\phi_m(\boldsymbol{r}-\boldsymbol{R})
 

\end{align*}
$$

と書くことができます。（一応、展開に使う関数にまだ制限は掛けていないので、まだ等式で書いています）

このように固体内の固有関数を、既知の関数系$\{\phi_m\}$で表現することができました。この展開を出発点として、展開係数を求めることができれば固有状態・固有エネルギーが分かります。

::: details （余談その2）
余談ですが、LCAO近似とTight-binding近似は割と混同されがちですが、別物であるようで、本章でもそれを意識した書き方にしています。
これは多分、原子軌道の線形結合を取ったとしても、十分空間的に広がった（束縛が弱い）電子の波動関数を表すこともできるし、また逆に強く束縛された電子の関数を表現するための方法は必ずしも原子軌道の線形結合でなくても良い、という事なのだと思いますが、詳細はわかりません。
とはいえ、ここで書いたことはこれから述べていく近似を追っていくことで具体的にイメージできるかと思います。
:::

## 固有値方程式の導出

さて、Bloch関数を原子軌道を用いて以下のように展開しました。

$$
\varphi_{\boldsymbol{k}}(\boldsymbol{r})

=
 \sum_{\boldsymbol{R}}e^{i\boldsymbol{k}\cdot\boldsymbol{R}} \sum_mb_m\phi_m(\boldsymbol{r}-\boldsymbol{R})
 
$$

Bloch関数は固体内の1体ハミルトニアンに対するシュレーディンガー方程式

$$
\hat{H}\varphi_{\boldsymbol{k}}(\boldsymbol{r}) = \varepsilon_{\boldsymbol{k}}\varphi_{\boldsymbol{k}}(\boldsymbol{r}),\\

\hat{H} = 
    -\frac{\hbar^2}{2m}\nabla^2
    +
   \sum_{\boldsymbol{R}}V(\boldsymbol{r} - \boldsymbol{R}),\\

\varphi_{\boldsymbol{k}}=\sum_{\boldsymbol{R}}e^{i\boldsymbol{k}\cdot\boldsymbol{R}} \sum_mb_m\phi_m(\boldsymbol{r}-\boldsymbol{R})
$$

を満たすので、あとはこの方程式を満たすための係数$b_m$を求めればBloch関数の原子軌道で展開した関数系を求めることができます。
ここでハミルトニアンが孤立原子のハミルトニアンではなく、孤立原子のポテンシャルを周期的に並べた固体のハミルトニアンであることに注意してください（念のため）

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
\varphi_{\boldsymbol{k}}(\boldsymbol{r})&=
\sum_mb_m
\left[\sum_{\boldsymbol{R}}e^{i\boldsymbol{k}\cdot\boldsymbol{R}} \phi_m(\boldsymbol{r}-\boldsymbol{R})
\right]

&\equiv
\sum_mb_m\Phi_{m,\boldsymbol{k}}(\boldsymbol{r})

\end{align*}
$$

と置き、上記シュレーディンガー方程式に左からBloch和$\Phi_n^*$を掛けて全空間について積分し、

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

を定義すると、最終的に固有エネルギー、そして係数$b_m$を求めるための方程式は、係数$b_m$を並べたベクトルを$\boldsymbol{b}$と置くと、以下の方程式

$$

  M_{\boldsymbol{k}}\boldsymbol{b} = \varepsilon_{\boldsymbol{k}} S_{\boldsymbol{k}} \boldsymbol{b}
$$

または行列式


$$
\begin{vmatrix}
  M_{\boldsymbol{k}} - \varepsilon_{\boldsymbol{k}} S_{\boldsymbol{k}}
\end{vmatrix} = 0
$$

と書くことができます。

なお、このように固有値が対角項以外にも入っているような方程式を、**一般化固有値方程式**や**一般化固有値問題**と呼ぶようです。行列$S_{\boldsymbol{k}}$が単位行列になれば、通常の固有値問題になります。また、上記の一般化固有値問題の解法は色々とあるようなので詳細が知りたい方はググってみてください。


さて、この方程式を解くことで（無限個の）固有値（エネルギーバンド）$\varepsilon_{n,\boldsymbol{k}}$と、対応した展開係数$b_m^n$、そしてその係数で展開されるBloch関数

$$
\varphi_{n,\boldsymbol{k}}(\boldsymbol{r}) = \sum_m b_m^n \sum_{\boldsymbol{R}}e^{i\boldsymbol{k}\cdot \boldsymbol{R}}\phi_m(\boldsymbol{r} - \boldsymbol{R})
$$

を得ることができます。



・・・とはいえ、これではほぼ、ただの行列表示のシュレーディンガー方程式で具体的にどのような原子軌道を展開に用いればよいのかが分から無ければ計算もできません。

というわけで次節でいよいよ、Tight-bindingな描像に基づいた近似を導入していき、現実的に計算可能な方程式を導いていきます。前置きが長すぎ？


# 第二工程：Tight-binding近似（前半）

長かった前置きが終わり、ここから「電子が格子点の局所ポテンシャルにタイトにバインディングされている」描像のもとで以下のような近似を導入し、計算可能な形で方程式を求めていきます。


## Tight-bindingな近似その1：バンド$n$の展開に必要な原子軌道の見積もり^[この部分はアシュクラフト・マーミン上Iを参考にした。]



まず初めに、重要な関係性として**各バンドに属する固有関数は少数の原子軌道関数で良く近似できる**ことを確かめます。^[私はこの部分が最近までイマイチ理解できておらず、Tight-bindingモデルの理解にかなり遠回りをしました]


さて、前節でWannier関数を原子軌道関数で（又はBloch関数をBloch和で）展開して、固有値・展開係数を決定するための一般化固有値問題展開を得ました。しかしながら「全ての原子軌道」を用意することは不可能なわけで、何らかの有限個の原子軌道のみの展開で上手く近似したいわけです。

そこで展開係数の関係を求めるために、まずは前節で得た固有値問題が首尾よく解けて、バンド$n$に属するBloch関数$\varphi_{n,\boldsymbol{k}}(\boldsymbol{r})$が得られたとしましょう。

この時そのBloch関数は


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

を満たします。ここに、今度は（Bloch和ではなく）原点を中心とした原子軌道関数の複素共役$\phi_{l}^*(\boldsymbol{r})$（$\boldsymbol{R} = \boldsymbol{0}$の場合に対応する）をかけて全空間について積分します。

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
d\boldsymbol{r}
\end{align*}
$$

特に$\varphi$が$\hat{H}$の固有状態のとき、固有値を$\varepsilon$とすると固有値は実数なので、$\int \left(\hat{H}\varphi(\boldsymbol{r})\right)^* \psi(\boldsymbol{r})d\boldsymbol{r} = \varepsilon\int \varphi^*(\boldsymbol{r}) \psi(\boldsymbol{r})d\boldsymbol{r}$、を利用して積分を計算していきます。

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

さて、先ほどの等式から$b_l^n$単体の項とそれ以外をそれぞれ移項したりして整理すると、展開係数$b_l^n$に関する等式

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
(\varepsilon_{n,\boldsymbol{k}} - \varepsilon_l^{\rm a})b_l^n 
$$

は小さくならなければならないことがわかります。これを満たすには

$$
\varepsilon_{n,\boldsymbol{k}} \simeq \varepsilon_l^{\rm a}
$$

すなわち固有値$\varepsilon_{n,\boldsymbol{k}}$に対して、係数$b_l^n$に対応した原子軌道関数の固有エネルギー$\varepsilon_l^{\rm a}$が近い値を持つ場合、または

$$
b_l^n \simeq 0
$$

の場合であると言えます。

これは、バンド$n$の固有エネルギーに近い原子準位を持つ原子軌道の展開係数は大きく、遠い原子準位の原子軌道の展開係数は小さいことを意味します。

従って、バンド$n$のBloch関数$\varphi_{n,\boldsymbol{k}}(\boldsymbol{r})$の展開は、そのBloch関数が属するエネルギーバンドの値と近い原子順位$\varepsilon_m^{\rm a}$を持つ（有限個の）原子軌道関数$\phi_m(\boldsymbol{r})$のみで良く近似される、ということを意味しています！

さらに、[固体（結晶）中の電子状態](https://zenn.dev/ponzumai/articles/tight-binding-model-electrons-in-solids)の章で見たように、バンド指標$n$に対応する固体内の電子の固有値は、
- バンド指標$n$ごとに離散的な、すなわち程度離れたエネルギーを持ち
- エネルギー$\varepsilon_{n,\boldsymbol{k}}$は連続、かつ限定的な領域「ブリルアンゾーン」を周期として持つ周期関数
- つまり上記固有値$\varepsilon_{n,\boldsymbol{k}}$は、BZ内のBloch波数\boldsymbol{k}$の関数として最大値と最小値を持ち、一定の幅に収まっている

わけでした。ここから、BZ内の全てのBloch波数$\boldsymbol{k}$に対して、「固有エネルギー$\varepsilon_{n,\boldsymbol{k}}$に近い原子準位」が共通で定義できると言えます。すなわち、$\varphi_{n,\boldsymbol{k}}(\boldsymbol{r})$はBloch波数$\boldsymbol{k}$によらず、共通の原子軌道のセットで展開できるわけです！！^[逆に、この近似の下ではある原子順位の重ね合わせで表現されるBloch関数のバンドは、概ねその原子準位に近い範囲で値を持つとも言えます。]


例えば縮退した、又は近い原子準位を持つ原子軌道関数を$m = m_1 , m_2,  m_3$として、そのグループをひとまず$q$でラベルしておき、原子準位を$\varepsilon^{{\rm a}_q}_{m_i}$と置くと、原子準位$\varepsilon^{{\rm a}_q}_{m_i}$に近い固有値に対応するBloch関数$\varphi_{n_q,\boldsymbol{k}}$は、3つの原子軌道（Bloch和）を用いて

$$
\varphi_{n_q,\boldsymbol{k}}(\boldsymbol{r}) \simeq \sum_{m = m_1,m_2,m_3}b_m^{n_q}
\sum_{\boldsymbol{R}}e^{i\boldsymbol{k}\cdot\boldsymbol{R}} \phi_m(\boldsymbol{r}-\boldsymbol{R})
$$

と近似できることを意味します。これは大きな進歩です。


### 余談その3
~~私のように~~余計なことが気になってしまうような人は、前掲の関係式

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
   \right]b_m^n \\

&\>\>\>\>+
\sum_m 
\left[
    \int \phi_l^*(\boldsymbol{r})
   \left(
   \sum_{\boldsymbol{R}'\neq \boldsymbol{0}}
      V(\boldsymbol{r}-\boldsymbol{R}')
      \right) \phi_m(\boldsymbol{r})d\boldsymbol{r}
      \right] b_m^n\\

&\>\>\>\>+
\sum_m 
\left[
   \sum_{\boldsymbol{R}\neq\boldsymbol{0}}
   e^{i\boldsymbol{k}\cdot\boldsymbol{R}} \int \phi_l^*(\boldsymbol{r})
   \left(
   \sum_{\boldsymbol{R}'\neq \boldsymbol{0}}
      V(\boldsymbol{r}-\boldsymbol{R}')
      \right)  \phi_m(\boldsymbol{r}-\boldsymbol{R})d\boldsymbol{r}
      \right]b_m^n
\end{align*}
$$

第1項について、「とはいえ、エネルギー順位が高い原子軌道関数は空間的にも広がっているだろうから、そういう場合を想定すれば第1項は大きくなり、「右辺が全部小さいから左辺も小さい、従ってエネルギーバンド近傍の原子順位の軌道だけで上手く近似できる」という結論にはならないのではないか？」等と~~いらんことを~~考えてしまうかもしれません。

そんな高エネルギー状態を展開に含めるのは物理的にはおかしいし、そもそも「局在した」関数で近似しよう、というのが出発点だとしたら上のようなことが気になるのもおかしいのかもしれませんが、一応安心するために、以下のような怪しい議論をしてみることにします。
**なお、以下のようなことを書いてある文献は私は寡聞にして知らず、まるで出鱈目なことを書いている可能性があるので注意してください**

（この後の章^[追記予定]で、さらに怪しい議論を展開して重なり積分が小さいことを納得していく予定です）

まず、式の簡略化のため、右辺第1項の「準位$l$の原子軌道関数$\phi_l$の重なり積分（と係数の積）の総和」を
$$
S_l^n\equiv\sum_m 
\left[
   \sum_{\boldsymbol{R}\neq\boldsymbol{0}}
   e^{i\boldsymbol{k}\cdot\boldsymbol{R}} 
   \int \phi_l^*(\boldsymbol{r})
   \phi_m(\boldsymbol{r}-\boldsymbol{R})d\boldsymbol{r}
   \right]b_m^n
$$

第2項、第3項を合わせた「原点以外の周期ポテンシャルを間に挟んだ積分（と係数の積）の相和」を

$$
V_l^n\equiv

\sum_m 
\left[
   \sum_{\boldsymbol{R}}
   e^{i\boldsymbol{k}\cdot\boldsymbol{R}} \int \phi_l^*(\boldsymbol{r})
   \left(
   \sum_{\boldsymbol{R}'\neq \boldsymbol{0}}
      V(\boldsymbol{r}-\boldsymbol{R}')
      \right)  \phi_m(\boldsymbol{r}-\boldsymbol{R})d\boldsymbol{r}
      \right]b^n_m 
$$

と置きます。

そして、バンド$n$と近い原子準位を持ち、したがって局在した（重なり積分が小さい）準位を$l$、バンドと離れ、大きい原子準位を持ち、したがって広がりの大きい（重なり積分の大きい）準位を$L$と置いて、展開係数$b_l$と$b_L$の比を求めてみます。


$$
\begin{align*}
\frac{b^n_L}{b^n_l} &= 
\frac{-S_L^n + (\varepsilon_{n,k} - \varepsilon_L^{\rm a})^{-1}V_L^n}{-S_l^n + (\varepsilon_{n,k} - \varepsilon_l^{\rm a})^{-1}V_l^n}\\

\end{align*}
$$

ここで、準位$l$の固有関数は局在性が強いと仮定したので、$S_l^n\simeq 0$とすると、

$$
\begin{align*}
\frac{b^n_L}{b^n_l} &\simeq
\frac{-S_L^n + (\varepsilon_{n,k} - \varepsilon_L^{\rm a})^{-1}V_L^n}{(\varepsilon_{n,k} - \varepsilon_l^{\rm a})^{-1}V_l^n}\\

&=
(\varepsilon_{n,k} - \varepsilon_l^{\rm a})\frac{-S_L^n }{V_l^n}
+\frac{\varepsilon_{n,k} - \varepsilon_l^{\rm a}}{\varepsilon_{n,k} - \varepsilon_L^{\rm a}} V_L^n

\end{align*}
$$

となります。ここで第2項も$\varepsilon_{n,k} \simeq \varepsilon_l^{\rm a}, \varepsilon_{n,k} \ll \varepsilon_L^{\rm a}$の仮定から小さくなり、無視して良さそうなので落として、最終的に「広がった」原子軌道と「局在した」原子軌道の展開係数の比

$$
\begin{align*}
\frac{b^n_L}{b^n_l} &\simeq

(\varepsilon_{n,k} - \varepsilon_l^{\rm a})\frac{-S_L^n }{V_l^n}

\end{align*}
$$

を得ます。この式を眺めてみると、例え重なり$S_L^n$が大きくても、周期ポテンシャルを含む積分部分


$$
V_l^n=

\sum_m 
\left[
   \sum_{\boldsymbol{R}}
   e^{i\boldsymbol{k}\cdot\boldsymbol{R}} \int \phi_l^*(\boldsymbol{r})
   \left(
   \sum_{\boldsymbol{R}'\neq \boldsymbol{0}}
      V(\boldsymbol{r}-\boldsymbol{R}')
      \right)  \phi_m(\boldsymbol{r}-\boldsymbol{R})d\boldsymbol{r}
      \right]b^n_m 
$$

がそれなりに大きければ、

$$
\begin{align*}
\frac{b^n_L}{b^n_l} &\simeq

(\varepsilon_{n,k} - \varepsilon_l^{\rm a})\frac{-S_L^n }{V_l^n}
\simeq 0

\end{align*}
$$

とでき、先ほどの結論に問題はなさそうです。
これは局所ポテンシャル（孤立原子のポテンシャル）が大きく、空間的に離れた原子軌道関数の減衰を補う位の積分値を与えてくれるような場合、という描像に対応してそうです。
「固体内の電子が格子点に強く束縛されている＝孤立原子のポテンシャルが強い」という、当初の「Tight-binding」な仮定とも一致しており、良い感じです。
（先ほどはこの部分を「小さい」と考えて同様の結論を導いたので矛盾しているように思えますが、比べている対象が違っていて、先ほどはエネルギー差$\varepsilon_{n,\boldsymbol{k}} - \varepsilon_l^{\rm a}$と比べて「小さい」としていて、今回は（広がった原子軌道の）重なり積分の相和$S_L^n$と比べていますが、展開係数$b_l^n$についての規格化条件を考えると、$S_L^n$は最大でも$1$ですので、それに比べて（それなりに大きい）という条件はそこまで矛盾しないかと思います。

さて、逆に周期ポテンシャルが弱い場合、極端な場合では重なり積分の和$S_l^n$とほぼ一致する場合は、$S_L^n/V_l^n$の分子が小さくなり、$(\varepsilon_{n,\boldsymbol{k}} - \varepsilon_l^{\rm a})$と打ち消しあっちゃって$b_L^n/b_l^n$が大きくなり得て、近似が上手くいかなさそうだと考えられます。

というわけで、それなりに強い局所ポテンシャルを考えている限りにおいては、安心してあるバンド$n$のBloch関数を少数のBloch和で（Wannier関数を少数の原子軌道関数で）展開できそうです。（まあそもそも、そこまでややこしいこと考えずに、展開してみて実験と合えばそれが正義なのかもしれませんが）

## 変分法による固有値方程式の（再）導出

以上のようにWannier関数が少数の原子軌道で、つまりBloch関数が少数のBloch和で展開されると近似した場合、Bloch関数$\varphi_{\boldsymbol{k}}(\boldsymbol{r})$の展開はエネルギー的に近い原子準位の集まりを用いてグループ分けできることがわかりました。

イメージ的には水素原子の固有関数のラベルを用いて、縮退した原子軌道を1グループとすると、

- $1s$グループ：$\varphi_{1s,\boldsymbol{k}}\simeq\sum_{\boldsymbol{R}}e^{i\boldsymbol{k}\cdot\boldsymbol{R}} \phi_{1s}(\boldsymbol{r}-\boldsymbol{R})$
- $2s$グループ：$\varphi_{2s,\boldsymbol{k}}\simeq\sum_{\boldsymbol{R}}e^{i\boldsymbol{k}\cdot\boldsymbol{R}} \phi_{2s}(\boldsymbol{r}-\boldsymbol{R})$
- $2p$グループ：$\varphi_{2p,\boldsymbol{k}}\simeq\sum_{m = 1,2,3}b_m^{2p}\sum_{\boldsymbol{R}}e^{i\boldsymbol{k}\cdot\boldsymbol{R}} \phi_{2p}(\boldsymbol{r}-\boldsymbol{R})$
- $3s$グループ：$\cdots$


みたいな感じ。~~（なお、この段階では異なる軌道間の重なり積分がまだ値を持つので、異なるグループに属する、つまり異なるバンドに属する（近似）ブロッホ関数は直交しません。後に重なり積分をゼロと置く近似をすることで異なるバンド間のBloch関数の正規直交性が満たされることになります）~~
↓
:::details （2023年1月20日追記）やっぱり直交するかも
エルミート演算子の異なる固有値を持つ固有関数は直交しなければならないので、Bloch和の重なり積分が値を持っていたとしても、近似した後のBloch関数は直交します（多分）。
これは図形的なベクトルでイメージすると、ある行列の固有ベクトルを求めたとして、その固有ベクトル同士は直交します。で、その固有ベクトルを、非直交の基底（例えば逆格子ベクトルを導入したときに考えた「斜め向きの」ベクトルみたいな）ベクトルで展開することもできるはずです。「固有ベクトルを展開した斜め向きの基底」は、（定義から明らかに）正規直交していませんが、だからと言ってそのベクトルで展開された固有ベクトルは、当然直交しているはずです。上の分で「固有ベクトル」を「固有関数」に、「非直交基底」を「Bloch和」に置き換えたような状況が、「異なるBloch和同士の重なり積分がゼロではない場合の固有関数の展開」の状況と一致するものだと思います。
ただ、全く別のBloch和で展開された固有関数同士の内積がどのようにしてゼロとなるのか（それぞれのブロック対角化された固有値方程式を解けば自動的に満たされるのか？）は、ちょっとまだよくわかってないので、いつかちゃんと考えて追記します。
:::

:::details （2023年1月20日追記その2）：やっぱり直交しなくてもいいかも
有限個のBloch和で近似した関数は、ハミルトニアンの正確な固有関数ではないので、べつに直交してなくてもいいような気もしてきました。
「エネルギー的に離れたBloch和同士だけ重なり積分を無視する」という近似に対応しているような気もしてきました。
頭が混乱しているので一度放置します。
:::

先述の通り、このように電子の固有状態を（少数の）原子軌道で展開する近似を**LCAO近似**と呼びます。

それではこのように展開できたとして、係数$b_m^n$はどのように決定すれば良いでしょうか？というと、それは結局、固体の1体ハミルトニアンに対するシュレーディンガー方程式に帰着されます。
それには一旦証明は省きますが、[He原子の固有状態を求めた際](https://zenn.dev/ponzumai/articles/tight-binding-model-many-electron-atom)にも用いた変分法を利用します。

すなわち、例えば$m = m_1, m_2, m_3$の原子準位に対応する原子軌道（このグループを先ほどと同様に$q$と名付けます）で展開した場合、つまり$\varepsilon_{n_q,\boldsymbol{k}} \simeq \varepsilon_{m_1}^{\rm a}, \varepsilon_{m_2}^{\rm a}, \varepsilon_{m_3}^{\rm a}$であるとき、Bloch関数

$$
\varphi_{n_q,\boldsymbol{k}}\simeq
\sum_{m}b_{m}^{n_q}
\left[\sum_{\boldsymbol{R}}e^{i\boldsymbol{k}\cdot\boldsymbol{R}} \phi^q_m(\boldsymbol{r}-\boldsymbol{R})
\right]
$$

を試行関数とし、「エネルギー期待値を最小にする」という条件を考えると、その条件は最終的には系のハミルトニアンの下でのシュレーディンガー方程式

$$
\hat{H}\varphi_{n_q,\boldsymbol{k}}(\boldsymbol{r}) = \varepsilon_{n_q,\boldsymbol{k}}\varphi_{n_q,\boldsymbol{k}}(\boldsymbol{r}),\\

\hat{H} = 
    -\frac{\hbar^2}{2m}\nabla^2
    +
   \sum_{\boldsymbol{R}}V(\boldsymbol{r} - \boldsymbol{R}),\\

\varphi_{n_q,\boldsymbol{k}}=
\sum_{m = m_1, m_2, m_3}b_m^{n_q}
\left[\sum_{\boldsymbol{R}}e^{i\boldsymbol{k}\cdot\boldsymbol{R}} \phi_m(\boldsymbol{r}-\boldsymbol{R})
\right]
$$

と一致することになります。また、最小エネルギーに限らず、変分法によって得られた2番目、3番目の固有値$\varepsilon_2^{\rm v}, \varepsilon_3^{\rm v}\cdots$（変分法"variation"の"v"）も、真の固有値$\varepsilon_2, \varepsilon_3,\cdots$と$\varepsilon_2^{\rm v} \leq \varepsilon_2,\>\>\>\>\varepsilon_3^{\rm v} \leq \varepsilon_3,\cdots$を満たすので（これも証明は一旦省きます）、上記「シュレーディンガー方程式」をから得られた複数の固有値は、そのまま系の近似的な固有値として考えることができます。

従ってこのグループの関数についてのBloch和を

$$
\Phi_{m,\boldsymbol{k}}^q\equiv \sum_{\boldsymbol{R}}e^{i\boldsymbol{k}\cdot\boldsymbol{R}}\phi_m(\boldsymbol{r} - \boldsymbol{R}), \>\>\>\>m = m_1, m_2, m_3
$$

として先ほどと同じくBloch和を掛けて行列

$$
\left(M_{\boldsymbol{k}}^q\right)_{ij} \equiv \int \Phi_{m_i,\boldsymbol{k}}^{q*}(\boldsymbol{r})\hat{H}\Phi^{q}_{m_j,\boldsymbol{k}}(\boldsymbol{r})d\boldsymbol{r},\\
\left(S_{\boldsymbol{k}}^q\right)_{ij} \equiv \int \Phi_{n,\boldsymbol{k}}^*(\boldsymbol{r})\Phi_{m,\boldsymbol{k}}(\boldsymbol{r})d\boldsymbol{r}
$$

を定義すると、最終的に固有エネルギー、そして係数$b_m^q$を求めるための方程式は以下の一般化固有値方程式で表されます。

$$
\begin{vmatrix}
  M_{\boldsymbol{k}}^q - \varepsilon_{\boldsymbol{k}} S_{\boldsymbol{k}}^q
\end{vmatrix} = 0
$$

この方程式を解くことで、展開に用いた関数の数だけ固有値・固有ベクトルが得られます。


上記のように、求めたいバンドに応じてエネルギー的に近そうな原子準位を選べば、その数程度のサイズの方程式（一般化固有値方程式）を解けばよいことになります。
とはいえ近似なので、完全に縮退した原子準位のみではなく、近いエネルギー順位も取り入れたりする必要はあるみたいですが、それにしても高々$4$個（$s+p$）とか、$6$個（$s+d$）位なので、十分計算できる範囲ですね。

あるいは、考えたい固体の構成要素が孤立原子の場合の内殻電子・価電子は知られているでしょうから、必要に応じて価電子やその近くの原子軌道のバンドを考えれば知りたい固体のエネルギーバンドを求めることができそうです。

さらに望むなら内殻電子から価電子まですべて考えたとしても、多くても十数～数十サイズの（一般化）固有値方程式にとどまるので、そのような固有値問題を数値的に解けば、十数～数十個のバンドを得ることができます。

これだけでもかなりの成果ですが、もう少し近似を考えていくことで、具体的な関数として固有エネルギーを得ることができます。


ただ、長くなってしまったのでここで一旦切り、続きは次章に回すことにします。


# おわりに

随分長くなってしまいましたが、本章では

- Bloch関数を原子軌道関数（Bloch和）で展開し、その際の固有値方程式$| M_{\boldsymbol{k}} - \varepsilon_{\boldsymbol{k}} S_{\boldsymbol{k}}| = 0$を求め
- あるバンドに属するBloch関数はお互いに近い原子準位を持つ少数の原子軌道関数のみを用いて近似的に展開$\varphi_{n_q,\boldsymbol{k}}(\boldsymbol{r}) \simeq \sum_{m = m_1,m_2,m_3}b_m^{n_q}\sum_{\boldsymbol{R}}e^{i\boldsymbol{k}\cdot\boldsymbol{R}} \phi_m(\boldsymbol{r}-\boldsymbol{R})$できることを示し
-最終的に上記少数の原子軌道で展開した場合の方程式を導出しました。


$$
\begin{vmatrix}
  M_{\boldsymbol{k}}^q - \varepsilon_{\boldsymbol{k}} S_{\boldsymbol{k}}^q
\end{vmatrix} = 0
$$


$$
\left(M_{\boldsymbol{k}}^q\right)_{ij} \equiv \int \Phi_{m_i,\boldsymbol{k}}^{q*}(\boldsymbol{r})\hat{H}\Phi^{q}_{m_j,\boldsymbol{k}}(\boldsymbol{r})d\boldsymbol{r},\\
\left(S_{\boldsymbol{k}}^q\right)_{ij} \equiv \int \Phi_{n,\boldsymbol{k}}^*(\boldsymbol{r})\Phi_{m,\boldsymbol{k}}(\boldsymbol{r})d\boldsymbol{r}
$$

これにより、例えば知りたい固体の孤立原子の内殻電子・価電子が分かっているような場合、価電子の（やそれと近い原子準位の）原子軌道のみでBlochを展開し、展開数に応じた（小さいサイズの）固有値問題を解くことで、価電子のエネルギーバンドを求めることができます。

後編では、ここからもう少しいくつかの近似を考え、具体的な固有エネルギー（エネルギーバンド）の関数系を求めていきます。