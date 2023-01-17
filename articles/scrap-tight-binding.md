---
title: "没原稿置き場"
emoji: "👏"
type: "tech" # tech: 技術記事 / idea: アイデア
topics: []
published: false
---
それには束縛条件の下での最適化法「ラグランジュの未定乗数法」を用いて、束縛条件は

$$
\begin{align*}
&\int\Phi^*(\tau_1,\tau_2)\Phi(\tau_1,\tau_2) d\tau_1d\tau_2\\
&= \int \varphi^*(1)\varphi^*(2)\varphi(1)\varphi(2)d\boldsymbol{r}_1d\boldsymbol{r}_2\\
&=\int\left|\varphi(\boldsymbol{r})\right|^2d\boldsymbol{r}\int\left|\varphi(\boldsymbol{r})\right|^2d\boldsymbol{r}\\
&=1
\end{align*}
$$

となり結局$\int\left|\varphi(\boldsymbol{r})\right|^2d\boldsymbol{r}=1$を満たせばよいので、未定乗数$\varepsilon$を用いて、次の式を最小にする関数$\varphi$を見つければよいことになります。

$$
\int \varphi^*(1)\varphi^*(2)\mathcal{H}\varphi(1)\varphi(2)d\boldsymbol{r}_1d\boldsymbol{r}_2 - \varepsilon
\left[
      \int\left|\varphi(\boldsymbol{r})\right|^2d\boldsymbol{r}-1  
\right]


$$
------------------------
$$
\begin{align*}
&=
 2\int
    \delta\varphi^*(\boldsymbol{r})
    \hat{H}
    \varphi(\boldsymbol{r})
    d\boldsymbol{r}
    +
    2\int
        \varphi^*(\boldsymbol{r})
    \hat{H}
    \delta\varphi(\boldsymbol{r})
    d\boldsymbol{r}
\end{align*}
$$



これは水素様原子中の電子のシュレーディンガー方程式

$$
\left(-\frac{\hbar^2}{2m} \nabla^2  -\frac{Ze^2}{4\pi\epsilon _0r}  \right) \varphi(\boldsymbol{r}) =\epsilon \varphi (\boldsymbol{r} )
$$

において第2項に$V(r)$が付け加わった形になっています。

したがって水素原子の場合で行ったように極座標表示に変換して変数分離をすることで、3つの微分方程式

$$
-\frac{\hbar^2}{2m} \left(
        \frac{d ^2R}{d r^2} + \frac{2}{r}\frac{d R}{d r}
        -\frac{\lambda }{r^2} R
        \right)
        +
        \left(
        -\frac{Ze^2}{4\pi\epsilon _0r}
        + V(r)
        \right)R = \epsilon R
       ,\\
\frac{1}{\sin\theta } \frac{d}{d\theta } \left( \sin\theta \frac{d\Theta }{d\theta }   \right) + \left( \lambda - \frac{m^2}{\sin^2\theta }  \right) \Theta =0,\\
\frac{d^2\Phi }{d\phi ^2} + m^2\Phi =0.
$$

と、クーロンポテンシャル部分に新たに平均場ポテンシャル$V(r)$が付け加わった形に整理できます。

ところで今回、一体の固有関数は球対称と仮定したのでした。そこで角度部分は定数、具体的には
$$
Y_0^0 = \frac{1}{\sqrt{4\pi}}
$$

となり、
試行関数$\varphi(|\boldsymbol{r}|)$を求めるための微分方程式は、動径部分のみの1変数微分方程式で$\lambda = 0$とした

面心立方格子の結晶内のポテンシャル

面心立方格子は以下の基本併進ベクトルで表され、

$$
\boldsymbol{a}_1 = \frac{a}{2}\boldsymbol{e}_x + \frac{a}{2}\boldsymbol{e}_y, \\ 
{}\\
\boldsymbol{a}_2 =  \frac{a}{2}\boldsymbol{e}_y + \frac{a}{2}\boldsymbol{e}_z,\\
{}\\
\boldsymbol{a}_3 =  \frac{a}{2}\boldsymbol{e}_x + \frac{a}{2}\boldsymbol{e}_z.
$$

結晶内のポテンシャルは格子ベクトルを

$$
\boldsymbol{n} = n_1\boldsymbol{a}_1 
+
n_2\boldsymbol{a}_2
+
n_3\boldsymbol{a}_3
$$

として$V(\boldsymbol{r}+\boldsymbol{n}) = V(\boldsymbol{r})$の周期を持ちます。この時逆格子ベクトルを

$$
\boldsymbol{K} = m_1\boldsymbol{b}_1 + m_2\boldsymbol{b}_2 + m_3\boldsymbol{b}_3,\\

\boldsymbol{b}_1 
=
\frac{\boldsymbol{a}_2\times\boldsymbol{a}_3}{\boldsymbol{a}_1\cdot(\boldsymbol{a}_2\times\boldsymbol{a}_3)}
=
\frac{2\pi}{a}(\boldsymbol{e}_x +\boldsymbol{e}_y - \boldsymbol{e}_z ) 
$$

と置くと

結晶内のポテンシャルは


$$
V(\boldsymbol{r}) = \sum_{\boldsymbol{K}}V_{\boldsymbol{K}}e^{i\boldsymbol{K}\cdot\boldsymbol{r}}\\
V_{\boldsymbol{K}}= \frac{1}{v_c}\int_{v_c} V(\boldsymbol{r})e^{-i\boldsymbol{K}\cdot\boldsymbol{r}}d\boldsymbol{r}
$$


と書けます。


$H_q$を対角化する行列$T$：

$$
T^{-1}H_qT = 
\begin{bmatrix}
\ddots &   &  &  &   \\
 & \>\varepsilon_{n,q} & \\
  &&\ddots  & & &\\
 & & & \>\>\varepsilon_{n',q} \\
 & &  &  & \ddots &  &  \\
   & & &&& \varepsilon_{n'',q} & \\
&& & & &&\ddots &  \\
\end{bmatrix}
$$


を考えると（ラベルのセンスが悪い気がしますが））、「固有ベクトル」

$$
T
\begin{bmatrix}
\vdots \\
c_{q - K_l} \\
\vdots \\
c_q \\
\vdots \\
c_{q + K_l} \\
\vdots
\end{bmatrix}

=
\begin{bmatrix}
\vdots \\
\sum_{m}T_{nm}c_{q - K_m} \\
\vdots \\
\sum_{m}T_{n'm}c_{q - K_m} \\
\vdots \\
\sum_{m}T_{n''m}c_{q - K_m} \\
\vdots
\end{bmatrix}
$$

と、（無限個の）$H_{q}$に対する$n$でラベルされた固有値

$$
\varepsilon_{n,q}
$$

を得ることができます。

（行列表示の都合で$K$を$K_m$、全ての$K$に関する和を$\sum_m$と書いていたところを、2行目で「全ての$K$に関する和$\sum_K$と書き換えました）


波数空間における周期性とワニエ関数
これは次の章に持っていくか



あるいは、もう少し物理的に考えてみると、
いわゆる「空格子近似」

逆格子ベクトルだけ異なる固有エネルギーは、波数の刻みに対して離散的だとみなせるので

ーーー

元々我々は多電子の固有状態を求めたいのでした。
そこでスレーター行列式を考えなおすと、固有状態
を「詰めて」行けば多電子の固有状態を得られるわけです。

今考えている周期的固有条件に対応して、N個の格子が存在していると考えられます。
格子一つ当り1個の電子がいると考えると、ブリルアンゾーン内にはN個の固有状態が並んでおり、その一つ一つの点に対して一つの1電子固有状態が対応します。

このような考えのもとエネルギー固有値が低い順からN個、1電子固有状態を詰めていったら、最高エネルギー固有値の値を考えることができます。そのようなエネルギーをフェルミエネルギーと呼び、$\varepsilon_F$と書きます。
この値はエネルギーバンド図では等エネルギー線として横軸と平行な直線で表されます。

バンド理論超概略

スレーター行列式と一粒子固有状態
フェルミエネルギー
エネルギーギャップ
　NFEモデルー摂動展開




Wannier関数参考文献：
https://cfm.ehu.es/ivo/publications/marzari-psik03.pdf
https://www.diva-portal.org/smash/get/diva2:631027/FULLTEXT01.pdf
http://www.ms.osakafu-u.ac.jp/~taguchi/densi/5th.pdf


# 第一量子化の最初の書き方：固有値方程式の理解がちょっとおかしかったのでいったんこちらに退避



# 第一工程：周期ポテンシャルの近似＋LCAO法による展開
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
w_{n,\boldsymbol{R}}(\boldsymbol{r}) = w_n(\boldsymbol{r}-\boldsymbol{R}) = \sum_mb_m^n\phi_m(\boldsymbol{r}-\boldsymbol{R}) 
$$

と展開します。（原子軌道関数であることがわかりやすいように、以降、原子軌道関数は先ほどと同じく$\phi$を用いて表します。）ここで$b_m^n$は規格化条件$\sum_m|b_m^n|^2 = 1$を満たすものとします。

このような考え方は、「原子軌道の線形結合」の英語バージョン"Linear Combination of Atomic Orbitals"の頭文字をとってLCAO法、またはLCAO近似（展開する原子軌道に制限をつけなければ「近似」ではないですが、多くの場合物理的な考察に応じて少数の原子軌道で展開することになります。（そうしないとそもそも計算できないし）そのように少数の原子軌道で展開するときに「LCAO近似」というのだと思います）と呼ばれています。

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

この係数を決定するために、過去にも用いた変分法を使います。すなわち、エネルギー期待値を最小化するような展開係数を求めていきます。

[He原子の固有状態を求めた](https://zenn.dev/ponzumai/articles/tight-binding-model-many-electron-atom)際には、試行関数として2電子のスレーター行列式を用いたことから変分法の結果平均場を含む方程式が出てきましたが、今回のような1電子のハミルトニアンの下では、「（試行）関数$\varphi_{n,\boldsymbol{k}}(\boldsymbol{r})$がエネルギー期待値を最小にする」という条件式は結局、元々考えていた固体内の1体ハミルトニアンに対するシュレーディンガー方程式

$$
\hat{H}\varphi_{n,\boldsymbol{k}}(\boldsymbol{r}) = \varepsilon_{n,\boldsymbol{k}}\varphi_{n,\boldsymbol{k}}(\boldsymbol{r}),\\

\hat{H} = 
    -\frac{\hbar^2}{2m}\nabla^2
    +
   \sum_{\boldsymbol{R}}V(\boldsymbol{r} - \boldsymbol{R}),\\

\varphi_{n,\boldsymbol{k}}=\sum_{\boldsymbol{R}}e^{i\boldsymbol{k}\cdot\boldsymbol{R}} \sum_mb_m^n\phi_m(\boldsymbol{r}-\boldsymbol{R})
$$

と一致します。ここでハミルトニアンが孤立原子のハミルトニアンではなく、孤立原子のポテンシャルを周期的に並べた固体のハミルトニアンであることに注意してください（念のため）

式を全部具体的に書くとこんな感じになります：

$$
\left(
 -\frac{\hbar^2}{2m}\nabla^2
    +
   \sum_{\boldsymbol{R}}V(\boldsymbol{r} - \boldsymbol{R})
   \right)
   \sum_{\boldsymbol{R}}e^{i\boldsymbol{k}\cdot\boldsymbol{R}} \sum_mb_m^n\phi_m(\boldsymbol{r}-\boldsymbol{R})
    = \varepsilon_{n,\boldsymbol{k}}\sum_{\boldsymbol{R}}e^{i\boldsymbol{k}\cdot\boldsymbol{R}} \sum_mb_m^n\phi_m(\boldsymbol{r}-\boldsymbol{R})
$$

ここで、Bloch関数の展開$\varphi_{n,\boldsymbol{k}}=\sum_{\boldsymbol{R}}e^{i\boldsymbol{k}\cdot\boldsymbol{R}} \sum_mb_m^n\phi_m(\boldsymbol{r}-\boldsymbol{R})$で和の順序を入れ替えてみると、

$$
\varphi_{n,\boldsymbol{k}}=
\sum_mb_m^n
\left[\sum_{\boldsymbol{R}}e^{i\boldsymbol{k}\cdot\boldsymbol{R}} \phi_m(\boldsymbol{r}-\boldsymbol{R})
\right]
$$

となり、これはBloch関数を、展開係数$b_m^n$、関数$\sum_{\boldsymbol{R}}e^{i\boldsymbol{k}\cdot\boldsymbol{R}} \phi_m(\boldsymbol{r}-\boldsymbol{R})$の線形結合で展開していると見ることもできます。なおこの関数はBloch sum（ブロッホ和、Bloch和）と呼ばれるみたいで、Blochの定理を満たします。また、一般に異なる格子点を中心とした原子軌道関数は直交しないので、ブロッホ和も直交しません。^[正規直交性を満たすようなブロッホ和の構成方法もあるみたいなのですが、この章では必要にならないので今はスルーしておきます。]

式を簡略化するためにBloch和を

$$
\begin{align*}
\varphi_{n,\boldsymbol{k}}&=
\sum_mb_m^n
\left[\sum_{\boldsymbol{R}}e^{i\boldsymbol{k}\cdot\boldsymbol{R}} \phi_m(\boldsymbol{r}-\boldsymbol{R})
\right]

&\equiv
\sum_mb_m^n\Phi_{m,\boldsymbol{k}}(\boldsymbol{r})

\end{align*}
$$

と置き、上記「シュレーディンガー方程式」に左からBloch和$\Phi_n^*$を掛けて全空間について積分し、

$$
\sum_mb_m^n\int \Phi_{n,\boldsymbol{k}}^*(\boldsymbol{r})\hat{H}\Phi_{m,\boldsymbol{k}}(\boldsymbol{r})d\boldsymbol{r}
=
\varepsilon_{n,\boldsymbol{k}}\sum_mb_m^n\int \Phi_{n,\boldsymbol{k}}^*(\boldsymbol{r})\Phi_{m,\boldsymbol{k}}(\boldsymbol{r})d\boldsymbol{r}
$$

行列

$$
\left(M_{\boldsymbol{k}}\right)_{nm} \equiv \int \Phi_{n,\boldsymbol{k}}^*(\boldsymbol{r})\hat{H}\Phi_{m,\boldsymbol{k}}(\boldsymbol{r})d\boldsymbol{r},\\
\left(S_{\boldsymbol{k}}\right)_{nm} \equiv \int \Phi_{n,\boldsymbol{k}}^*(\boldsymbol{r})\Phi_{m,\boldsymbol{k}}(\boldsymbol{r})d\boldsymbol{r}
$$

を定義すると、最終的に固有エネルギー、そして係数$b_m^n$を求めるための方程式は以下の固有値方程式で表され、Bloch波数$\boldsymbol{k}$をパラメータとして$n$バンド分の固有値・固有ベクトルが得られます。

$$
\begin{Vmatrix}
  M_{\boldsymbol{k}} - \varepsilon_{\boldsymbol{k}} S_{\boldsymbol{k}}
\end{Vmatrix} = 0
$$

・・・とはいえ、これではただの行列表示のシュレーディンガー方程式で、このままで解けるなら何の苦労もないわけで、Wannier関数（Bloch関数）の展開に選んだ関数の数だけある行列の対角化が必要になっちゃいます。

というわけで次節でいよいよ、Tight-bindingな描像に基づいた近似を導入していきます。前置きが長すぎ？

# 第二工程：Tight-binding近似

長かった前置きが終わり、ここから「電子が格子点の局所ポテンシャルにタイトにバインディングされている」描像のもとで以下のような近似を導入し、具体的な固有エネルギー、固有関数を求めていきます。

- 

## Tight-bindingな近似その1：バンド$n$の展開に必要な原子軌道の見積もり^[この部分はアシュクラフト・マーミン上（I)を参考にした。]^[ものすごく余談ですが、最初にゼミで読んだ翻訳書がニジェ・オーランドでGrassmann代数をさらっとかじった後アルトランド・サイモンズで経路積分をかじってさらに例題を解いてくみたいな内容だったので、カナカナ連名著者の教科書は自分には敷居の高い高尚な内容が書かれているものだと思い込んでアシュクラフト・マーミン等も避けていたのですが、最近になって読んでみると非常に丁寧に知りたいことが書いてあり愕然としました。さっさと読めばよかった。。。]



まず初めに、重要な関係性として**各バンド$n$に属する固有関数は少数の原子軌道関数で良く近似できる**ことを確かめます。^[私はこの部分が最近までイマイチ理解できておらず、Tight-bindingモデルの理解にかなり遠回りをしました]


さて、前節でWannier関数を原子軌道関数で（又はBloch関数をBloch和で）展開したのですが、展開に用いる原子軌道が多くなればなるほど精度が良くなるように思えます。それはそうなのですが、そんなことを言うと近似も近似でなくなるわけで、以下のようにある係数$b_m$、つまり展開に含まれる原子軌道のうち、$m$番目のエネルギー順位の原子軌道が含まれる割合についての関係式を得ることができ、結局そのエネルギーバンド程度の固有エネルギーに対応する（少数の）原子軌道関数のみで展開できることが分かります。

早速前に進みましょう。行列形式は一旦忘れて、先ほど変分法によって初めに得られた、固有状態が満たすべき方程式（＝シュレーディンガー方程式）

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

に、今度は原点を中心とした原子軌道関数の複素共役$\phi_{l}^*(\boldsymbol{r})$（$\boldsymbol{R} = \boldsymbol{0}$の場合に対応する）をかけて全空間について積分します。

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
なお、ここで現れる積分

$$
\int \phi_l^*(\boldsymbol{r})
   \phi_m(\boldsymbol{r}-\boldsymbol{R})d\boldsymbol{r}
$$

は**Overlap Integral（重なり積分）**と呼ばれます。

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

なのでした。つまりバンド指標$n$で指定されるの固有エネルギーは、BZ内のBloch波数$\boldsymbol{k}$に応じて値に幅を持ちますが、最大値と最小値を持ち、一定の幅に収まっているわけです。

従って、「あるバンド$n$のBloch関数$\varphi_{n,\boldsymbol{k}}(\boldsymbol{r})$またはWannier関数$w_{n,\boldsymbol{R}}(\boldsymbol{r})$の展開に用いる原子軌道関数$\phi_l(\boldsymbol{r})$」は、そのBloch関数が属するエネルギーバンドの幅と近い原子順位$\varepsilon_m^{\rm a}$を持つ原子軌道関数$\phi_m$のみで良く近似される、ということを意味しています！

逆に、この近似の下ではある原子順位の重ね合わせで表現されるBloch関数のバンドは、概ねその原子準位に近い範囲で値を持つと言えます。

もう少し具体的に書くと、バンド$n$に属するBloch関数$\varphi_{n,\boldsymbol{k}}$は、縮退した、又は近い原子準位を持つ原子軌道関数を$m = m_1 , m_2,  m_3$とすると、

$$
\varphi_{n,\boldsymbol{k}}(\boldsymbol{r}) \simeq \sum_{m = m_1,m_2,m_3}b_m^n
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

第1項について、「とはいえ、エネルギー順位が高い原子軌道関数は空間的にも広がっているだろうから、そういう場合を想定すれば第1項は大きくなり、「右辺が全部小さいから左辺も小さい、従ってエネルギーバンド近傍の原子順位の軌道だけで上手く近似できる」という結論にはならないのではないか？」等と~~いらんことを~~考えてしまうかもしれません。

そんな高エネルギー状態を展開に含めるのは物理的にはおかしいし、そもそも「局在した」関数で近似しよう、というのが出発点だとしたら上のようなことが気になるのもおかしいのかもしれませんが、一応安心するために、以下のような怪しい議論をしてみることにします。
**なお、以下のようなことを書いてある文献は私は寡聞にして知らず、まるで出鱈目なことを書いている可能性があるので注意してください**


まず、式の簡略化のため、右辺第1項の「準位$l$の原子軌道関数$\phi_l$の重なり積分（と係数の積）の相和」を
$$
S_l^n\equiv\sum_m 
\left[
   \sum_{\boldsymbol{R}\neq\boldsymbol{0}}
   e^{i\boldsymbol{k}\cdot\boldsymbol{R}} 
   \int \phi_l^*(\boldsymbol{r})
   \phi_m(\boldsymbol{r}-\boldsymbol{R})d\boldsymbol{r}
   \right]b^n_m 
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

そして、バンド$n$と近い原子準位を持ち、したがって局在した（重なり積分が小さい）準位を$l$、バンドと離れ、大きい原子準位を持ち、したがって広がりの大きい（重なり積分の大きい）準位を$L$と置いて、展開係数$b^n_l$と$b^n_L$の比を求めてみます。


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

となります。ここで第2項も$\varepsilon_{n,k} \simeq \varepsilon_l^{\rm a}, \varepsilon_{n,k} \ll \varepsilon_L^{\rm a}$の仮定から無視して良さそうなので落として、最終的に「広がった」原子準位と「局在した」原子準位の比

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

さて、逆に周期ポテンシャルが弱い場合、極端な場合では重なり積分の和$S_l^n$とほぼ一致する場合は、$S_L^n/V_l^n$の分子が小さくなり、$(\varepsilon_{n,k} - \varepsilon_l^{\rm a})$と打ち消しあっちゃって$b_L^n/b_l^n$が大きくなり得て、近似が上手くいかなさそうだと考えられます。

というわけで、それなりに強い局所ポテンシャルを考えている限りにおいては、安心してバンド$n$のBloch関数を少数のBloch和で（Wannier関数を少数の原子軌道関数で）展開できそうです。（まあそもそも、そこまでややこしいこと考えずに、展開してみて実験と合えばそれが正義なのかもしれませんが）

### 固有値方程式のブロック対角化

以上のようにWannier関数が少数の原子軌道で、つまりBloch関数が少数のBloch和で展開されると仮定した場合、Bloch関数$\varphi_{n,\boldsymbol{k}}(\boldsymbol{r})$はバンド指標$n$について、サイズの小さい行列でブロック対角化された形になります。すなわち、バンド$n$について例えば$m = m_1, m_2, m_3$の原子準位に対応する原子軌道で展開される場合、つまり$\varepsilon_{n,\boldsymbol{k}} \simeq \varepsilon_{m_1}^{\rm a}, \varepsilon_{m_2}^{\rm a}, \varepsilon_{m_3}^{\rm a}$であるとき、Bloch関数の展開は

$$
\varphi_{n,\boldsymbol{k}}=
\sum_{m = m_1, m_2, m_3}b_m^n
\left[\sum_{\boldsymbol{R}}e^{i\boldsymbol{k}\cdot\boldsymbol{R}} \phi_m(\boldsymbol{r}-\boldsymbol{R})
\right]
$$

となるため、工程1で求めた固有値方程式は対応する原子準位$\phi_{m_i}(\boldsymbol{r}), i = 1,2,3$を用いて、次の$3\times 3$行列の行列式

$$
\begin{vmatrix}
(M_{\boldsymbol{k}})_{11} -\varepsilon_{n,\boldsymbol{k}}(S_{\boldsymbol{k}})_{11} & (M_{\boldsymbol{k}})_{12} -\varepsilon_{n,\boldsymbol{k}}(S_{\boldsymbol{k}})_{12} & (M_{\boldsymbol{k}})_{13} -\varepsilon_{n,\boldsymbol{k}}(S_{\boldsymbol{k}})_{13} \\

(M_{\boldsymbol{k}})_{21} -\varepsilon_{n,\boldsymbol{k}}(S_{\boldsymbol{k}})_{21} & (M_{\boldsymbol{k}})_{22} -\varepsilon_{n,\boldsymbol{k}}(S_{\boldsymbol{k}})_{22} & (M_{\boldsymbol{k}})_{23} -\varepsilon_{n,\boldsymbol{k}}(S_{\boldsymbol{k}})_{23} \\

(M_{\boldsymbol{k}})_{31} -\varepsilon_{n,\boldsymbol{k}}(S_{\boldsymbol{k}})_{31} & (M_{\boldsymbol{k}})_{32} -\varepsilon_{n,\boldsymbol{k}}(S_{\boldsymbol{k}})_{32} & (M_{\boldsymbol{k}})_{33} -\varepsilon_{n,\boldsymbol{k}}(S_{\boldsymbol{k}})_{33} \\
\end{vmatrix}
=0
$$

を解けばよいことになります。凄い省エネです。


＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝

