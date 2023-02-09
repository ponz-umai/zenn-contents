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





例えば$m_1,m_2,m_3$の3つの原子軌道で展開したとすると、

$$
\frac{M_{\boldsymbol{k}}}{N} =

\begin{bmatrix}
\varepsilon_{m_1}^{\rm a} - \Delta\varepsilon_{m_1m_1} & \sum_{\boldsymbol{R}_I}e^{i\boldsymbol{k}\cdot\boldsymbol{R}_I} t_{\boldsymbol{R}_I}^{m_1m_2} -\Delta\varepsilon_{m_1m_2} & \sum_{\boldsymbol{R}_I}e^{i\boldsymbol{k}\cdot\boldsymbol{R}_I} t_{\boldsymbol{R}_I}^{m_1m_3}- \Delta\varepsilon_{m_1m_3}\\

\sum_{\boldsymbol{R}_I}e^{i\boldsymbol{k}\cdot\boldsymbol{R}_I} t_{\boldsymbol{R}_I}^{m_2m_1} - \Delta\varepsilon_{m_2m_1} &
  \varepsilon_{m_2}^{\rm a} -\Delta\varepsilon_{m_2m_2} &
   \sum_{\boldsymbol{R}_I}e^{i\boldsymbol{k}\cdot\boldsymbol{R}_I} t_{\boldsymbol{R}_I}^{m_2m_3} - \Delta\varepsilon_{m_2m_3}\\

\sum_{\boldsymbol{R}_I}e^{i\boldsymbol{k}\cdot\boldsymbol{R}_I} t_{\boldsymbol{R}_I}^{m_3m_1} - \Delta\varepsilon_{m_3m_1} &
 \sum_{\boldsymbol{R}_I}e^{i\boldsymbol{k}\cdot\boldsymbol{R}_I} t_{\boldsymbol{R}_I}^{m_3m_2} - \Delta\varepsilon_{m_3m_2} & 
 \varepsilon_{m_3}^{\rm a} - \Delta\varepsilon_{m_3m_3}
\end{bmatrix}
$$


$$
\\

\begin{bmatrix}
0 & \sum_{\boldsymbol{R}_I}e^{i\boldsymbol{k}\cdot\boldsymbol{R}_I} t_{\boldsymbol{R}_I}^{m_1m_2}  &
 \sum_{\boldsymbol{R}_I}e^{i\boldsymbol{k}\cdot\boldsymbol{R}_I} t_{\boldsymbol{R}_I}^{m_1m_3}\\

\sum_{\boldsymbol{R}_I}e^{i\boldsymbol{k}\cdot\boldsymbol{R}_I} t_{\boldsymbol{R}_I}^{m_2m_1} &
  0 &
   \sum_{\boldsymbol{R}_I}e^{i\boldsymbol{k}\cdot\boldsymbol{R}_I} t_{\boldsymbol{R}_I}^{m_2m_3} \\

\sum_{\boldsymbol{R}_I}e^{i\boldsymbol{k}\cdot\boldsymbol{R}_I} t_{\boldsymbol{R}_I}^{m_3m_1}  &
 \sum_{\boldsymbol{R}_I}e^{i\boldsymbol{k}\cdot\boldsymbol{R}_I} t_{\boldsymbol{R}_I}^{m_3m_2}  & 
 0
\end{bmatrix}\\

\Rightarrow

\begin{vmatrix}
\varepsilon & \sum_{\boldsymbol{R}_I}e^{i\boldsymbol{k}\cdot\boldsymbol{R}_I} t_{\boldsymbol{R}_I}^{m_1m_2}  &
 \sum_{\boldsymbol{R}_I}e^{i\boldsymbol{k}\cdot\boldsymbol{R}_I} t_{\boldsymbol{R}_I}^{m_1m_3}\\

\sum_{\boldsymbol{R}_I}e^{i\boldsymbol{k}\cdot\boldsymbol{R}_I} t_{\boldsymbol{R}_I}^{m_2m_1} &
  \varepsilon &
   \sum_{\boldsymbol{R}_I}e^{i\boldsymbol{k}\cdot\boldsymbol{R}_I} t_{\boldsymbol{R}_I}^{m_2m_3} \\

\sum_{\boldsymbol{R}_I}e^{i\boldsymbol{k}\cdot\boldsymbol{R}_I} t_{\boldsymbol{R}_I}^{m_3m_1}  &
 \sum_{\boldsymbol{R}_I}e^{i\boldsymbol{k}\cdot\boldsymbol{R}_I} t_{\boldsymbol{R}_I}^{m_3m_2}  & 
 \varepsilon
\end{vmatrix}=0,\\

\varepsilon \equiv \varepsilon_{\boldsymbol{k}} - (\varepsilon^{\rm a} - \Delta\varepsilon)

$$


^[今、$s$軌道は実数関数、つまり$\phi_s(\boldsymbol{r})^* = \phi_s(\boldsymbol{r})$で、かつ角度部分は定数となり$r$のみに依存します。$V$は引力ポテンシャルであることから$V(\boldsymbol{r})<0$です。原子軌道が局在している仮定から、積分の値の大部分は、原点から$\boldsymbol{R}$離れた部分の寄与が大きくなるはずです。ここで動径部分が、原点と、$\boldsymbol{r} = \boldsymbol{R}$で符号を変えないのであれば、積分全体は負になりそうですが、別にそうとも限らない気がするので、$s$バンドの重なり積分が必ず負になるかどうかは私はまだよくわかっていません。]




# 飛び移り積分関連

$$
\hat{H}^{\rm c} \varphi_{n,\boldsymbol{k}}(\boldsymbol{r}) = \varepsilon_{n,\boldsymbol{k}}\varphi_{n,\boldsymbol{k}}(\boldsymbol{r})
$$

は、孤立原子ハミルトニアン$\hat{H}^{\rm a} = \frac{-\hbar^2}{2m}\nabla^2 + V(\boldsymbol{r})$に対して、固有値（原子準位）$\varepsilon_m^{\rm a}$を持つ固有状態である原子軌道を$\phi_m(\boldsymbol{r})$からなるBloch和$\Phi_m(\boldsymbol{r}) = \sum_{\boldsymbol{R}}e^{i\boldsymbol{k}\cdot\boldsymbol{R}} \phi_m(\boldsymbol{r} - \boldsymbol{R})$のうち、固体電子の固有値と原子準位が近い（$\varepsilon_{n,\boldsymbol{k}} \sim \varepsilon_m^{\rm a}$）軌道を集めて線形結合を取った形で以下のように近似できます（ここでは軌道の数を$M$個とします）：

$$
\varphi_{n,\boldsymbol{k}}(\boldsymbol{r}) \simeq
\sum_m b_m^n\Phi_m(\boldsymbol{r}) ,\\

\Phi_{m,\boldsymbol{k}}(\boldsymbol{r}) = \sum_{\boldsymbol{R}}e^{i\boldsymbol{k}\cdot\boldsymbol{R}} \phi_m(\boldsymbol{r} - \boldsymbol{R}).
$$

また、上の展開式を逆に考えると、

$$
(B)_{nm} = b_{m}^{n}
$$

から行列$B$を定義し、また$B^{-1}B = I$（$I$は単位行列）を満たす逆行列$B^{-1}$の要素を

$$
(B^{-1})_{nm} \equiv \tilde{b}_{m}^{n}
$$

と書くことにして、これらと$\phi_m(\boldsymbol{r}) = N^{-1}\sum_{\boldsymbol{k}}\Phi_m(\boldsymbol{r})$


または上の表現を$-\boldsymbol{R}$だけずらすと、飛び移り積分

$$
\begin{align*}
t_{\boldsymbol{R}_I}^{m'm} &= 
-
\int_V\phi_{m'}^*(\boldsymbol{r})
  V(\boldsymbol{r}-\boldsymbol{R}_I)
\phi_m(\boldsymbol{r}-\boldsymbol{R}_I)d \boldsymbol{r} \\

&=
-
\int_V\phi_{m'}^*(\boldsymbol{r} + \boldsymbol{R}_I)
  V(\boldsymbol{r})
\phi_m(\boldsymbol{r})d \boldsymbol{r} 
\end{align*}
$$

は（隣接）格子ベクトル$\boldsymbol{R}_I$で指定される格子点に局在した状態$m$の電子が、原点に、原子軌道$\phi_{m'}(\boldsymbol{r}+\boldsymbol{R})$となって飛び移ってくる確率、とも言えます。

# 飛び移り積分の物理的意味：飛び移り積分の定義変更前

展開係数$c_{m',\boldsymbol{R}}^m$は、$\phi_{m'}(\boldsymbol{r} - \boldsymbol{R})$の係数が$c_{m',-\boldsymbol{R}}^m$であることに注意して、

$$
c_{m',\boldsymbol{R}}^m = \int\phi^*_{m'}(\boldsymbol{r} + \boldsymbol{R})\left(
\frac{-\hbar^2}{2m}\nabla^2 + \sum_{\boldsymbol{R}'}V(\boldsymbol{r}-\boldsymbol{R}')
\right)\phi_m(\boldsymbol{r})d\boldsymbol{r}
$$

となります。さらに、孤立原子ハミルトニアン部分の固有関数であることを利用して、

$$
\begin{align*}
（上式右辺）&=

\varepsilon_m^{\rm a}\int\phi_{m'}^*(\boldsymbol{r} + \boldsymbol{R})\phi_m(\boldsymbol{r})d \boldsymbol{r} + 
\int\phi_{m'}^*(\boldsymbol{r} + \boldsymbol{R})
\left(
  \sum_{\boldsymbol{R}' \neq -\boldsymbol{R}}V(\boldsymbol{r} - \boldsymbol{R}')
\right)
\phi_m(\boldsymbol{r})d \boldsymbol{r} \\

&=
\varepsilon_m^{\rm a}\delta_{m,m'}\delta_{\boldsymbol{R},\boldsymbol{0}}+ 
\int\phi_{m'}^*(\boldsymbol{r} + \boldsymbol{R})
\left(
  \sum_{\boldsymbol{R}' \neq -\boldsymbol{R}}V(\boldsymbol{r} - \boldsymbol{R}')
\right)
\phi_m(\boldsymbol{r})d \boldsymbol{r} 

\end{align*}
$$

と書き換えることができます。ついに「飛び移り積分」っぽい形の積分がでてきましたね。

ここで$\boldsymbol{R} = \boldsymbol{0}$とその他で分けると、

### $\boldsymbol{R} = \boldsymbol{0}$

$$
\begin{align*}
c_{m',\boldsymbol{0}}^m &\simeq

\varepsilon_m^{\rm a}\delta_{mm'}

+
\int\phi_{m'}^*(\boldsymbol{r})
\left(
  \sum_{\boldsymbol{R}' \neq \boldsymbol{0}}V(\boldsymbol{r} - \boldsymbol{R}')
\right)
\phi_m(\boldsymbol{r})d \boldsymbol{r} \\

&=
\varepsilon_m^{\rm a}\delta_{mm'}

+\Delta\varepsilon_{m'm}
\end{align*}
$$

と、原子準位$\varepsilon_m^{\rm a}$と結晶場積分$\Delta\varepsilon_{m'm}$で書けることがわかります。

### $\boldsymbol{R} \neq \boldsymbol{0}$

また$\boldsymbol{R} \neq \boldsymbol{0}$の場合は、前章と同様にTight-binding近似として、（最）隣接格子間の飛び移り積分ではないものをゼロと置くと、


$$
c_{m',\boldsymbol{R}}^m \simeq

\delta_{\boldsymbol{R},\boldsymbol{R}_I}\sum_{\boldsymbol{R}_I} 
\int\phi_{m'}^*(\boldsymbol{r} + \boldsymbol{R}_I)
  V(\boldsymbol{r})
\phi_m(\boldsymbol{r})d \boldsymbol{r} 
$$

を得ます。ここで積分$\int\phi_{m'}^*(\boldsymbol{r} + \boldsymbol{R}_I)V(\boldsymbol{r})\phi_m(\boldsymbol{r})d \boldsymbol{r}$は、ちゃんと書けば積分範囲が結晶全体の定積分ですが、積分変数を$\boldsymbol{r} +\boldsymbol{R}_I\rightarrow \boldsymbol{r}$と置きなおすと、周期的境界条件より積分範囲は変わらず、飛び移り積分の形

$$
\int_V\phi_{m'}^*(\boldsymbol{r})
  V(\boldsymbol{r}-\boldsymbol{R}_I)
\phi_m(\boldsymbol{r}-\boldsymbol{R}_I)d \boldsymbol{r} \equiv - t_{\boldsymbol{R}_I}^{m'm}
$$

と一致します。以上より、展開係数は

$$
c_{m',\boldsymbol{R}}^m \simeq
(\varepsilon_m^{\rm a}\delta_{mm'} +\Delta\varepsilon_{m'm})\delta_{\boldsymbol{R},\boldsymbol{0}}

-
\delta_{\boldsymbol{R},\boldsymbol{R}_I} t_{\boldsymbol{R}_I}^{m'm}
$$

となり、原子軌道関数に結晶のハミルトニアンを作用させた関係式は

$$
\begin{align*}
\hat{H}^{\rm c} \phi_m(\boldsymbol{r})
&=
\left(
\frac{-\hbar^2}{2m}\nabla^2 + \sum_{\boldsymbol{R}}V(\boldsymbol{r}-\boldsymbol{R})
\right) \phi_m(\boldsymbol{r}) \\


&\simeq
\sum_{m,\boldsymbol{R}}
\left\{
  (\varepsilon_m^{\rm a}\delta_{mm'} +\Delta\varepsilon_{m'm})\delta_{\boldsymbol{R},\boldsymbol{0}}
  -
  \delta_{\boldsymbol{R},\boldsymbol{R}_I} t_{\boldsymbol{R}_I}^{m'm}
\right\}
\phi_{m}(\boldsymbol{r} - \boldsymbol{R})\\

&=
(\varepsilon_m^{\rm a} + \Delta\varepsilon_{mm})\phi_m(\boldsymbol{r})

+
\sum_{m'} \Delta\varepsilon_{m'm}\phi_{m'}(\boldsymbol{r})\\

&\>\>\>\>+\sum_{\boldsymbol{R}_I}\sum_{m'}(-t_{\boldsymbol{R}_I}^{m'm})\phi_{m'}(\boldsymbol{r} + \boldsymbol{R}_I)

\end{align*}
$$

となることが分かりました。（$\phi_{m'}(\boldsymbol{r} -\boldsymbol{R})$の係数を$c_{m'-\boldsymbol{R}}^m$と置いたので、$\boldsymbol{R}$の符号がプラスになっています。）


・・・だからどうしたという感じですが、ここで記憶の片隅から、波動関数にハミルトニアンを作用させた意味を引っ張り出してきます。



## 波動関数の時間発展

さて、定常状態ばっかり扱っていてもはや記憶の彼方に葬り去られていましたが、本来のシュレーディンガー方程式は時間と位置に関する微分方程式

$$
i\hbar \frac{\partial}{\partial t}\psi(\boldsymbol{r},t) = \hat{H}\psi(\boldsymbol{r},t)
$$

なのでした。さらに、時刻$t$から微小時間$\Delta t$の変化を考えると、左辺の偏微分を

$$
\frac{\partial}{\partial t}\psi(\boldsymbol{r},t)\simeq

\frac{\psi(\boldsymbol{r},t + \Delta t) - \psi(\boldsymbol{r},t)}{\Delta t}
$$

として、シュレーディンガー方程式の両辺を$i\hbar$で割って整理すると、微小時間についての波動関数の時間発展の式

$$
\frac{\psi(\boldsymbol{r},t + \Delta t) - \psi(\boldsymbol{r},t)}{\Delta t}

\simeq

\frac{-i}{\hbar}\hat{H}\psi(\boldsymbol{r},t)\\

\Rightarrow
\psi(\boldsymbol{r},t + \Delta t) \simeq 

\psi(\boldsymbol{r},t ) + \frac{-i}{\hbar}\Delta t\hat{H}\psi(\boldsymbol{r},t)
$$

を得ることができます。これで準備が整ったので、**結晶の中、$t=0$で原点に一つだけ孤立原子軌道がある場合の時間発展**を考えてみましょう。すなわちハミルトニアン、$t=0$の波動関数を

$$
\begin{align*}
&\hat{H} = \hat{H}^{\rm c} 
=
\left(
\frac{-\hbar^2}{2m}\nabla^2 + \sum_{\boldsymbol{R}}V(\boldsymbol{r}-\boldsymbol{R})
\right),\\
&\psi(\boldsymbol{r},0) = \phi_m(\boldsymbol{r})
\end{align*}
$$

の時間発展を考えます。これは先ほど得た関係式より、


$$
\begin{align*}
\psi(\boldsymbol{r},\Delta t) &\simeq 

\psi(\boldsymbol{r},0 ) + \frac{-i}{\hbar}\Delta t\hat{H}\psi(\boldsymbol{r},0)
\\

&\simeq
\phi_m(\boldsymbol{r})
+
\frac{-i}{\hbar}\Delta t\hat{H}\phi_m(\boldsymbol{r})
\\

&\simeq

\left\{1-\frac{-i}{\hbar}\Delta t
(\varepsilon_m^{\rm a} + \Delta\varepsilon_{mm})
\right\}
\phi_m(\boldsymbol{r})\\

&\>\>\>\>+
\frac{-i}{\hbar}\Delta t\sum_{m'} \Delta\varepsilon_{m'm}\phi_{m'}(\boldsymbol{r})\\

&\>\>\>\>+\frac{-i}{\hbar}\Delta t\sum_{\boldsymbol{R}_I}\sum_{m'}(-t_{\boldsymbol{R}_I}^{m'm})\phi_{m'}(\boldsymbol{r} + \boldsymbol{R}_I)


\end{align*}
$$

が得られます。すなわち、原点にポツンと1つ原子軌道があるとき、微小時間後の波動関数を各格子点に平行移動した、様々な準位の原子軌道$\{\phi_{m'}(\boldsymbol{r} + \boldsymbol{R}\}$の重ね合わせで書くことができました。この時の$\phi_{m'}(\boldsymbol{r}+\boldsymbol{R})$の係数はかなり雑に言うと「電子が格子点$-\boldsymbol{R}$上で、$m'$の原子軌道で存在する確率」つまり「原点の格子点$\boldsymbol{R} = \boldsymbol{0}$にいた状態$m$の電子が、格子点$-\boldsymbol{R}$に状態$m'$で「飛び移る」確率」と考えられそうです。

#### 注意点1
ここで位置を表す変数が$\boldsymbol{r}$と、$\boldsymbol{R}$の二つ出てきて紛らわしいのですが、前者$\boldsymbol{r}$は「電子の波動関数の形」$\simeq$「電子の確率雲の形」に対応する変数で、後者$\boldsymbol{R}$は「波動関数の中心」$\simeq$「固体の中の電子の場所」に対応する変数です。この表現もかなり正確性がアヤシイですが。。。

#### 注意点2

正確（？）には、波動関数は観測することができず、観測できるのはエルミート演算子で表される何かしらの物理量です。そして、ある波動関数を何らかの演算子の固有状態で展開していた時、物理量の測定を通してその固有状態である波動関数が確定します（なんかこの書き方もイマイチな気がする・・）。そこで、波動関数$\phi_m(\boldsymbol{r}-\boldsymbol{R})$は何か物理量に対応する演算子の固有状態になっているかというと、そこが良くわかっていないので、例えば上記のような時間発展を考えたとしても、時刻$\Delta t$後に、「電子が格子点$\boldsymbol{R}$上で、$m'$の原子軌道で存在する」かどうかを確認できるかというとそういうわけではないように思えます。

## 原子軌道の飛び移り積分の物理的意味


特に、微小時間$\Delta t$後に、電子が原点から（隣接）格子ベクトル$-\boldsymbol{R}_I$で指定される格子点に局在した原子軌道$\phi_{m'}(\boldsymbol{r}+\boldsymbol{R})$となっているような確率は、飛び移り積分$t_{\boldsymbol{R}_I}^{m'm}$：

$$
\begin{align*}
t_{\boldsymbol{R}_I}^{m'm} &= 
-
\int_V\phi_{m'}^*(\boldsymbol{r})
  V(\boldsymbol{r}-\boldsymbol{R}_I)
\phi_m(\boldsymbol{r}-\boldsymbol{R}_I)d \boldsymbol{r} \\

&=
-
\int_V\phi_{m'}^*(\boldsymbol{r} + \boldsymbol{R}_I)
  V(\boldsymbol{r})
\phi_m(\boldsymbol{r})d \boldsymbol{r} 
\end{align*}
$$

に比例することがわかります。ここで、最後の式では$\phi_{m'}(\boldsymbol{r}+\boldsymbol{R})$と座標をそろえるために再度積分変数を変更しました。



以上のように（きわめて怪しい）議論を繰り広げることで、飛び移り積分はその名前の通り、**原点（あるいはある格子点$\boldsymbol{R}'$）にいる$m$状態の原子軌道$\phi_m(\boldsymbol{r})$（$\phi_m(\boldsymbol{r}-\boldsymbol{R}')$）が、（微小時間後に）状態$m'$になって原点（あるいはある格子点$\boldsymbol{R}'$）から格子ベクトル$-\boldsymbol{R}_I$離れた場所へ、$\phi_{m'}(\boldsymbol{r}+\boldsymbol{R}_I)$（$\phi_{m'}(\boldsymbol{r}-\boldsymbol{R}'+\boldsymbol{R}_I)$）として飛び移って行く確率**（に比例する量）を表すことが分かりました。

ただ、飛び移り積分の（教科書通りの）定義を採用すると、「飛び移り先」が$-\boldsymbol{R}_I$になってしまうのが気になるところなのです。何か対称性をちゃんと考えれば解消されるか、どこかで勘違いしているかもしれませんが、一旦思考停止しておきます。この後来る第二量子化とかを考えていけば何か解消されるかもしれません。

~~まあ、これがわかったところで、微小時間後の波動関数が測定できるわけでもないし、何の意味があるのかは分かりませんが。というか自己満足ですが。~~


# 重なり積分






一応、前章の展開

$$
\begin{align*}
\hat{H}^{\rm c}\phi_m(\boldsymbol{r} - \boldsymbol{R})&=
\sum_{m'}
\sum_{\boldsymbol{R}'}

\left(\sum_{n}\varepsilon_{n,\boldsymbol{R}-\boldsymbol{R}'}\tilde{b}_m^n b_{m'}^n
\right)\phi_{m'}(\boldsymbol{r} - \boldsymbol{R}')


\\
&\equiv
\sum_{m',\boldsymbol{R}'}c_{m',\boldsymbol{R}'}^{m,\boldsymbol{R}}\phi_{m'}(\boldsymbol{r}-\boldsymbol{R}').
\end{align*}
$$

において、重なり積分をゼロとしなければ、展開係数

$$
\begin{align*}
c_{m',\boldsymbol{R}'}^{m,\boldsymbol{R}} 
&= \int \phi_{m'}(\boldsymbol{r} - \boldsymbol{R}')\hat{H}^{\rm c}\phi_m(\boldsymbol{r} - \boldsymbol{R})d\boldsymbol{r} 

- 
\sum_{m',\boldsymbol{R}'\neq \boldsymbol{R}}c_{m',\boldsymbol{R}'}^{m,\boldsymbol{R}}s_{\boldsymbol{R},\boldsymbol{R}'}^{m,m'}\\

&\simeq

-t_{(m',\boldsymbol{R}')\leftarrow (m,\boldsymbol{R})}

-
\sum_{m',\boldsymbol{R}'\neq \boldsymbol{R}}c_{m',\boldsymbol{R}'}^{m,\boldsymbol{R}}s_{\boldsymbol{R},\boldsymbol{R}'}^{m,m'}
\end{align*}
$$

のような感じで「飛び移り積分の補正」みたいな形で出てきますが、そもそも別の展開係数が掛かってしまってますし、この方向で考えると前章の内容が全部成立しなくなっちゃうので、上手い方向ではなさそうです。






エネルギー固有値に対応する固有状態：原点の原子軌道

たまたま隣の原子に束縛された状態の電子が、ふわふわと漂っていたとする

その電子が、原点の原子核ポテンシャルに束縛される確率

一見するとポテンシャルに束縛されてもいないし、飛び移り積分よりも大きな値（確率）を持ちそうである。

しかし飛び移り積分は、全方向を格子ポテンシャルに囲まれており、逃げ場はない。とどまるか、どこかに移るか。

一方原子核ポテンシャルは、ちょうどよい角度でいかないと、変な方向に飛んで行ったり、束縛しきれなかったりするだろう

というわけで小さいのも納得である



重なり積分
一見、「広がっている」から、「重なり」は大きい気もする

でも物理的な描像で考えてみると、「違う原子に束縛された瞬間に大きなエネルギーを持ってしまう確率

これは小さいと考えられるし、実際、広がっていたとしても、実際は重なりは小さいのだと納得できるような気もする



＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝











多電子の波動関数は、軌道部分の波動関数$\varphi_a, \varphi_b\cdots$と、上向き/下向きスピン状態$\varphi_{1\uparrow}, \varphi_{1\downarrow},\varphi_{2\uparrow}, \varphi_{2\downarrow}\cdots$、または別の表現をするとスピン関数$\alpha(\sigma)$、$\beta(\sigma)$の積を考え、軌道部分の関数をラベルする量子数（例えば水素様原子中の電子の場合だと$(n,l,m)$）とスピン状態のラベル$\uparrow,\downarrow$をまとめて、ギリシャ文字$\lambda,\mu,\xi$などと表した「スピン軌道関数」$\varphi_\lambda(\tau_1), \varphi_\mu(\tau_2),\cdots,\varphi_\xi(\tau_N)$を用いて、

$$
\begin{align*}
\Phi_\mathrm{F}(\tau_1, \tau_2,\cdots,\tau_N) &= 
\frac{1}{\sqrt{N!}}
\begin{vmatrix}
\varphi_\lambda(\tau_1) & \varphi_\mu(\tau_1) & \cdots & \varphi_\xi(\tau_1)\\
\varphi_\lambda(\tau_2) & \varphi_\mu(\tau_2) & \cdots & \varphi_\xi(\tau_2)\\
& \cdots & \cdots\\
\varphi_\lambda(\tau_N) & \varphi_\mu(\tau_N) & \cdots & \varphi_\xi(\tau_N)
\end{vmatrix}\\

&=\frac{1}{\sqrt{N!}}\sum_P(-1)^P\hat{P}\varphi_\lambda(\tau_1)\varphi_\mu(\tau_2)\cdots\varphi_\xi(\tau_{N})\\

&\equiv \left| \varphi_a \overline{\varphi}_b\cdots \overline{\varphi}_n\right|
\end{align*}
$$

と書けることを[この章](https://zenn.dev/ponzumai/articles/tight-binding-model-spin)で確認しました。左辺二つ目の式は、行列式を、座標を入れ替える置換演算子$\hat{P}$と、置換の回数だけ$(-1)$を掛ける$(-1)^P$で表したものです。また多体波動関数の添え字$\rm F$は、フェルミオンの$\rm F$です。

ここで、軌道部分の波動関数$\varphi_a, \varphi_b\cdots$について正規直交性を満たすものとし、また完全系をなすと仮定します。（スピン関数もスピン空間内で完全系を成すので、その積のスピン軌道関数も完全系をなすと仮定していることになります）

というわけで多体ハミルトニアンをスレーター行列式に作用させてみると、

$$
\begin{align*}
&\mathcal{H}\Phi_{\rm F}(\tau_1, \tau_2,\cdots,\tau_N)\\



&=
\left|
\left(\hat{H}\varphi_\lambda\right)\varphi_\mu\cdots\varphi_\xi
    \right|

+
\left|
\varphi_\lambda\left(\hat{H}\varphi_\mu\right)\cdots\varphi_\xi
    \right|

+\cdots

+
\left|
\varphi_\lambda\varphi_\mu\cdots\left(\hat{H}\varphi_\xi\right)
    \right|

\end{align*}
$$

と、スレーター行列式の各1電子の波動関数のどれか一つに対して、1体ハミルトニアンが作用したスレーター行列式の総和となります。

:::details 証明

証明というか確認ですが、実際にスレーター行列式に作用させた後に、出てきた項を並べ替えることで以下のように示せます。

$$
\begin{align*}
&\mathcal{H}\Phi_{\rm F}(\tau_1, \tau_2,\cdots,\tau_N)\\

&=
\sum_i\hat{H}_i
\frac{1}{\sqrt{N!}}\sum_P(-1)^P\hat{P}\varphi_\lambda(\tau_1)\varphi_\mu(\tau_2)\cdots\varphi_\xi(\tau_{N})\\

&=
\frac{1}{\sqrt{N!}}
\left[
\left\{
    \hat{H}_1+ \hat{H}_2 + \cdots\hat{H}_N
    \right\}
    \left\{
(-1)^0\varphi_\lambda(\tau_1)\varphi_\mu(\tau_2)\cdots\varphi_\xi(\tau_{N})
        \right. \right.\\
        
&\>\>\>\>\>\>\>\>\>\>\>\>\>\>\>\>\>\>\>\>

+(-1)^1\varphi_\lambda(\tau_2)\varphi_\mu(\tau_1)\cdots\varphi_\xi(\tau_{N})
+\cdots
        \left.\left.
        
        \right\}
    \right]\\

&=
\frac{1}{\sqrt{N!}}
\left[
    (-1)^0\left(\hat{H}_1\varphi_\lambda(\tau_1)\right)\varphi_\mu(\tau_2)\cdots\varphi_\xi(\tau_{N}) 
    +
    (-1)^0\varphi_\lambda(\tau_1)\left(\hat{H}_2\varphi_\mu(\tau_2)\right)\cdots\varphi_\xi(\tau_{N}) 
    +\cdots
    \right.
    \\
&\>\>\>\>\>\>\>\>\>\>\>\>\>\>\>\>\>
\left.
+
(-1)^1\varphi_\lambda(\tau_2)\left(\hat{H}_1\varphi_\mu(\tau_1)\right)\cdots\varphi_\xi(\tau_{N})
+
(-1)^1\left(\hat{H}_2\varphi_\lambda(\tau_2)\right)\varphi_\mu(\tau_1)\cdots\varphi_\xi(\tau_{N})
+\cdots
\right]\\

&=
\frac{1}{\sqrt{N!}}
\left[
    (-1)^0\left(\hat{H}_1\varphi_\lambda(\tau_1)\right)\varphi_\mu(\tau_2)\cdots\varphi_\xi(\tau_{N}) 
    +
    (-1)^1\left(\hat{H}_2\varphi_\lambda(\tau_2)\right)\varphi_\mu(\tau_1)\cdots\varphi_\xi(\tau_{N})
    + \cdots\right]
\\

&\>\>\>\>+
\frac{1}{\sqrt{N!}}
\left[

    (-1)^0\varphi_\lambda(\tau_1)\left(\hat{H}_2\varphi_\mu(\tau_2)\right)\cdots\varphi_\xi(\tau_{N}) 
    
+
(-1)^1\varphi_\lambda(\tau_2)\left(\hat{H}_1\varphi_\mu(\tau_1)\right)\cdots\varphi_\xi(\tau_{N})
+

+\cdots
\right]\\

&=
\left|
\left(\hat{H}\varphi_\lambda\right)\varphi_\mu\cdots\varphi_\xi
    \right|

+
\left|
\varphi_\lambda\left(\hat{H}\varphi_\mu\right)\cdots\varphi_\xi
    \right|

+\cdots

+
\left|
\varphi_\lambda\varphi_\mu\cdots\left(\hat{H}\varphi_\xi\right)
    \right|.

\end{align*}
$$
:::

# TB第二量子化表示






この場合格子定数を$a$とすると、基本格子ベクトルは


$$
\boldsymbol{a}_1 = (a,0,0),\\
\boldsymbol{a}_2 = (0,a,0),\\
\boldsymbol{a}_3 = (0,0,a)
$$

となり、隣接格子ベクトルは基本格子ベクトルの逆方向まで考えて

$$
\boldsymbol{N}_i = (\pm a,0,0), (0,\pm a,0),(0,0,\pm a)
$$

です。また、再隣接への飛び移り積分は、$x,y,z$方向が等価なので全ての再隣接ベクトルに対して同じ値になるため、

$$
\begin{align*}
-t_{(s,\boldsymbol{R}+\boldsymbol{N}_i) \leftarrow (s,\boldsymbol{R})}
&=

\int
 \phi_s(|\boldsymbol{r}-\boldsymbol{N}_i|)
 
   V(|\boldsymbol{r} - \boldsymbol{N}_i|)
 
 \phi_s(|\boldsymbol{r}|)d\boldsymbol{r}\\


&=

\int
 \phi_s(|\boldsymbol{r}-(\pm a,0,0)|)
 
   V(|\boldsymbol{r} - (\pm a,0,0)|)
 
 \phi_s(|\boldsymbol{r}|)d\boldsymbol{r}=\cdots\\

&\equiv

-t_s^{\rm NN}


\end{align*}
$$

と置けます。（"NN"は、Nearlest Neighbour"の頭文字）

このように考えるとハミルトニアンは飛び移り積分を和の外に出せて、



$$
\begin{align*}
\mathcal{H}^s

&=
        -t_{s}^{\rm NN} 
\sum_{\gamma=\uparrow,\downarrow}
    \sum_{\boldsymbol{R}}
    \sum_{i\in I}


    \hat{a}_{s,\boldsymbol{R}+\boldsymbol{N}_i,\gamma}^\dagger\hat{a}_{s,\boldsymbol{R},\gamma}
    
\end{align*}
$$

となります。これはかなり見覚えのある形になりました。ここで和$\sum_{\boldsymbol{R}}\sum_{i\in I}$ですが、



＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝
## s軌道の具体例





=========================================

ここで、第一量子化で考えたように、「固体の固有関数をとりあえず3次元格子

が、今、1つの（$s$軌道の）原子軌道関数のみを考えているということは、固体の固有関数をこの原子軌道関数のみのBloch和で近似できると考えていることに対応します。（あるいは（Wannier関数を一つの原子軌道関数で近似しているとも見れます）

この時固有関数を$\varphi_{s,\boldsymbol{k}}(\boldsymbol{r})$とすると


$$
\varphi_{s,\boldsymbol{k}}(\boldsymbol{r})\simeq
\sum_{\boldsymbol{R}}e^{i\boldsymbol{k}\cdot\boldsymbol{R}}\phi_s(\boldsymbol{r}-\boldsymbol{R})
$$

$$
\phi_s(\boldsymbol{r}-\boldsymbol{R})
=
\frac{1}{N^3}
\sum_{\boldsymbol{k}}e^{-i\boldsymbol{k}\cdot\boldsymbol{R}}
\varphi_{s,\boldsymbol{k}}(\boldsymbol{r})
$$

特に、$\boldsymbol{R} = n\boldsymbol{a}_1$の時は、

$$
\phi_s(\boldsymbol{r}-n\boldsymbol{a}_1)
=
\frac{1}{N^3}
\sum_{\boldsymbol{k}}
e^{-ik_1na}
\varphi_{s,\boldsymbol{k}}(\boldsymbol{r})
$$

これを変換と考えると
$$
\hat{a}^\dagger_{s,j,\gamma}

=
\sum_{\boldsymbol{k}} \braket{ \varphi_{s,\boldsymbol{k},\gamma}|\phi_{s,j,\gamma}}\hat{b}_{s,\boldsymbol{k},\gamma}^\dagger
=
\frac{1}{N^3}\sum_{\boldsymbol{k}}e^{-ik_1aj}b_{s,\boldsymbol{k},\gamma}^\dagger
$$

$$
\hat{a}_{s,j,\gamma} = \sum_{\boldsymbol{k}} \braket{\phi_{s,j,\gamma} | \varphi_{s,\boldsymbol{k},\gamma}}\hat{b}_{s,\boldsymbol{k},\gamma}

= \sum_{\boldsymbol{k}}e^{ik_1aj}b_{s,\boldsymbol{k},\gamma}
$$


となる。これを代入すると、

$$
\begin{align*}
\mathcal{H}^{s,1D}

&=
   -t_{s}^{\rm NN}
\sum_{\gamma=\uparrow,\downarrow}
    \sum_{j}
    \sum_{\delta=\pm 1}

  \frac{1}{N^3}\sum_{\boldsymbol{k},\boldsymbol{k}'}
  e^{i((k_1-k_1')a_1\times j)}
    e^{-ik_1'a_1\delta}
  b_{s,\boldsymbol{k}',\gamma}^\dagger
  b_{s,\boldsymbol{k},\gamma}
    \\

&=
   -t_{s}^{\rm NN} \frac{1}{N^3}
\sum_{\gamma=\uparrow,\downarrow}
    \sum_{\delta=\pm 1}

 \sum_{\boldsymbol{k},\boldsymbol{k}'}
  N\delta_{k_1,k_1'}
  e^{-ik_1'a_1\delta}
  b_{s,\boldsymbol{k}',\gamma}^\dagger
  b_{s,\boldsymbol{k},\gamma}
    \\

&=
   -t_{s}^{\rm NN} \frac{1}{N^2}
\sum_{\gamma=\uparrow,\downarrow}
 \sum_{k_1}\sum_{k_2,k_3,k_2',k_3'}
 \sum_{\delta=\pm 1}
  e^{-ik_1a_1\delta}
  b_{s,\boldsymbol{k},\gamma}^\dagger
  b_{s,\boldsymbol{k},\gamma}
    \\

&=
   -t_{s}^{\rm NN} \frac{1}{N^2}
\sum_{\gamma=\uparrow,\downarrow}
 \sum_{k_1}
 \sum_{k_2,k_3,k_2',k_3'}
 2\cos(k_1a_1)
  b_{s,\boldsymbol{k},\gamma}^\dagger
  b_{s,\boldsymbol{k},\gamma}
    \\
\end{align*}
$$

ここでスレーター行列式$|\varphi_{\boldsymbol{k}}\cdots|$に対して

$$
2\cos(k_1a_1)
\sum_{k_2,k_3,k_2',k_3'}
 
  b_{s,\boldsymbol{k}',\gamma}^\dagger
  b_{s,\boldsymbol{k},\gamma}
  |\varphi_{\boldsymbol{k}}\cdots|

=
2\cos(k_1a_1)
\sum_{k_2,k_3}
 
  b_{s,\boldsymbol{k},\gamma}^\dagger
  b_{s,\boldsymbol{k},\gamma}

  |\varphi_{\boldsymbol{k}}\cdots|
$$



＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝



$$
\hat{a}^\dagger_{s,n + \delta,\gamma}
=
\frac{1}{N^3}
\sum_{\boldsymbol{k}} e^{-i(k_1a_1\times (n + \delta))}
b_{s,\boldsymbol{k},\gamma}^\dagger,\\

\hat{a}_{s,j,\gamma} = 
\sum_{\boldsymbol{k}} 
e^{i(k_1a_1\times j)}
b_{s,\boldsymbol{k},\gamma}
$$

ここで

$$
\hat{a}_\nu = \sum_\eta c_{\eta,\nu}^* \hat{b}_\eta = \sum_\eta \braket{\phi_\eta|\varphi_\nu} \hat{b}_\eta
$$





==================================




$$
\varphi_{s,k_1}(\boldsymbol{r})\simeq
\sum_{j}e^{ik_1a\times j}\phi_s(\boldsymbol{r}-j\boldsymbol{a}_1)
$$

逆変換

$$
\phi_s(\boldsymbol{r}-j\boldsymbol{a}_1)
=
\frac{1}{N}
\sum_{k_1}
e^{-ik_1a\times j}
\varphi_{s,k_1}(\boldsymbol{r})
$$

であると考えられます。


これを基底$\phi_s(\boldsymbol{r}-j\boldsymbol{a}_1)$から$\varphi_{s,k_1}(\boldsymbol{r})$への変換と考え、$\varphi_{s,k_1}(\boldsymbol{r})$とスピン関数の積の生成消滅演算子を$\hat{b}_{s,k_1,\gamma}^\dagger,\hat{b}_{s,\boldsymbol{k},\gamma}$とすれば、生成消滅演算子の変換は前後でスピン状態は変わらないので

$$
\hat{a}^\dagger_{s,j + \delta,\gamma}
=
\frac{1}{N}
\sum_{k_1} e^{-i(k_1a_1\times (j + \delta))}
b_{s,k_1,\gamma}^\dagger,\\

\hat{a}_{s,j,\gamma} = 
\frac{1}{N}
\sum_{k_1} 
e^{i(k_1a_1\times j)}
b_{s,k_1,\gamma}
$$

となります。これをハミルトニアンに代入すると、


$$
\begin{align*}
\mathcal{H}^{s,1D}

&=
   -t_{s}^{\rm NN}
\sum_{\gamma=\uparrow,\downarrow}
    \sum_{j}
    \sum_{\delta=\pm 1}

  \frac{1}{N^6}\sum_{\boldsymbol{k},\boldsymbol{k}'}
  e^{i((k_1-k_1')a_1\times j)}
    e^{-ik_1'a_1\delta}
  b_{s,\boldsymbol{k}',\gamma}^\dagger
  b_{s,\boldsymbol{k},\gamma}
    \\

&=
   -t_{s}^{\rm NN} \frac{1}{N^6}
\sum_{\gamma=\uparrow,\downarrow}
    \sum_{\delta=\pm 1}

 \sum_{\boldsymbol{k},\boldsymbol{k}'}
  N\delta_{k_1,k_1'}
  e^{-ik_1'a_1\delta}
  b_{s,\boldsymbol{k}',\gamma}^\dagger
  b_{s,\boldsymbol{k},\gamma}
    \\

&=
   -t_{s}^{\rm NN} \frac{1}{N^5}
\sum_{\gamma=\uparrow,\downarrow}
 \sum_{k_1}\sum_{k_2,k_3,k_2',k_3'}
 \sum_{\delta=\pm 1}
  e^{-ik_1a_1\delta}
  b_{s,\boldsymbol{k},\gamma}^\dagger
  b_{s,\boldsymbol{k},\gamma}
    \\

&=
   -t_{s}^{\rm NN} \frac{1}{N^6}
\sum_{\gamma=\uparrow,\downarrow}
 \sum_{k_1}
 \sum_{k_2,k_3,k_2',k_3'}
 2\cos(k_1a_1)
  b_{s,\boldsymbol{k},\gamma}^\dagger
  b_{s,\boldsymbol{k},\gamma}
    \\
\end{align*}
$$



## 1次元の具体例


この時、大きな特徴として



ですが、$s$軌道関数が球対称実関数、局所ポテンシャル$V(\boldsymbol{r})$も球対称実関数と置いているので、どちらの方向の飛び移り積分も等しくかつ実数の定数、

$$
\begin{align*}
-t_{(s,\boldsymbol{R}\pm\boldsymbol{a}_1) \leftarrow (s,\boldsymbol{R})}
&=

\int
 \phi_s(|\boldsymbol{r}\pm\boldsymbol{a}_1|)
 
   V(|\boldsymbol{r} \pm \boldsymbol{a}_1|)
 
 \phi_s(|\boldsymbol{r}|)d\boldsymbol{r}\\

 &\equiv-t_s^{\rm NN} = -(t_s^{\rm NN})^*\\

 
\end{align*}
$$

となります。これらを踏まえて第二量子化表示のハミルトニアンは


$$
\begin{align*}
\mathcal{H}^{s,1D}

&=
   
\sum_{\gamma=\uparrow,\downarrow}
    \sum_{\boldsymbol{R}}
    \sum_{i\in I}

     -t_{(s,\boldsymbol{R}+\boldsymbol{N}_i) \leftarrow ({s},\boldsymbol{R})}
    \hat{a}_{s,\boldsymbol{R}+\boldsymbol{N}_i,\gamma}^\dagger\hat{a}_{s,\boldsymbol{R},\gamma}
    \\
&\simeq

\sum_{\gamma=\uparrow,\downarrow}
    \sum_{\boldsymbol{R}}
    \sum_{\boldsymbol{N}_i=\pm \boldsymbol{a}_1}

     -t_{s}^{\rm NN}
    \hat{a}_{s,\boldsymbol{R}+\boldsymbol{N}_i,\gamma}^\dagger\hat{a}_{s,\boldsymbol{R},\gamma}
    \\


\end{align*}
$$


となります。ここで、このハミルトニアンはスレーター行列式の中の1電子状態$\phi_{s,\boldsymbol{R}}$のラベル$\boldsymbol{R}$を、$\boldsymbol{a}_1$方向にしか変えないため、スレーター行列式を1次元上に並んだ原子軌道関数$\phi_s(\boldsymbol{r}), \phi_s(\boldsymbol{r}\pm \boldsymbol{a}_1),\phi_s(\boldsymbol{r}\pm 2\boldsymbol{a}_1),\cdots ,\phi_s(\boldsymbol{r}\pm i\boldsymbol{a}_1)\cdots$のみを考えることにして^[この辺の説明がまだ自分でも混乱している気がするのでいつか修正する]、関数

$$
\phi_s(\boldsymbol{r}\pm i\boldsymbol{a}_1)
$$

（とスピン関数の積）の生成消滅演算子を

$$
    \hat{a}_{s,i,\gamma}^\dagger,\hat{a}_{s,i,\gamma}
$$

と書くことにすると、1次元的なハミルトニアン

$$
\begin{align*}
\mathcal{H}^{s,1D}

&=
   -t_{s}^{\rm NN}
\sum_{\gamma=\uparrow,\downarrow}
    \sum_{i}
    \sum_{\delta=\pm 1}

  
    \hat{a}_{s,i+\delta,\gamma}^\dagger\hat{a}_{s,i,\gamma}
    \\


\end{align*}
$$

が得られます。ここで、$\sum_{i}\sum_{\delta=\pm 1}$ですが、全ての$i$とその隣接項$i\pm 1$についての和を取るので、これを隣接した$i,j$のみ和を取るという意味で、$i<j$として

$$
-t_{s}^{\rm NN}
    \sum_{<i,j>}
    \left(\hat{a}_{s,i,\gamma}^\dagger\hat{a}_{s,j,\gamma}
+
\hat{a}_{s,j,\gamma}^\dagger\hat{a}_{s,i,\gamma}
\right)
$$

と書くこともあります。
余談ですが、この第2項をエルミート共役の意味で"+h.c."と書いていることが多いですが、エルミート共役だと$\hat{a}_{s,j,\gamma}\hat{a}_{s,i,\gamma}^\dagger = -\hat{a}_{s,i,\gamma}^\dagger\hat{a}_{s,j,\gamma}$になってしまうのでつじつまが合わない気がするのですが、正直この辺はあまりよくわかっていません。

# 水素原子

動径関数の一般項。これ自分がどこ見て書いたかすら覚えていない

$$
R_{nl}(r)= e^{-r/2\alpha}F(\alpha r) ,\alpha^2 =\frac{8m|\epsilon |}{\hbar^2} , \\
F(\rho ) = \rho ^lL_{n+l}^{2l+1},\\
L_{n+l}^{2l+1} = \sum_{k=0}^{n-l-1}(-1)^{k+2l+1}\frac{\left\{ (n+l)! \right\}^2 }{(n-l-1-k)!(2l+1+k)!k!}\rho ^k ,\\
\Theta(\theta ) \Phi(\phi )= Y_l^m(\theta ,\phi ),\\

$$

# 飛び移り積分の物理的意味

## 複素共役について間違ったことを書いていた当り

さて、最後に複素共役の謎を解明しておきます。これまでは話の展開を簡単にするために、原点にいる原子軌道の時間発展を考えてきました。するとなぜか、飛び移り積分の複素共役が出てきてしまいました。

次に一般の、格子点$\boldsymbol{R}$に局在した原子軌道$\phi_m(\boldsymbol{r} - \boldsymbol{R})$の時間発展を考えてみますと、

$$
\hat{H}^{\rm c}\phi_m(\boldsymbol{r} - \boldsymbol{R})
$$

は、「証明」部分と同じような計算をすれば、原点に局在した原子軌道の場合と同様に


$$
\begin{align*}
\hat{H}^{\rm c}\phi_m(\boldsymbol{r} - \boldsymbol{R})&=
\sum_{m'}
\sum_{\boldsymbol{R}'}

\left(\sum_{n}\varepsilon_{n,\boldsymbol{R}-\boldsymbol{R}'}\tilde{b}_m^n b_{m'}^n
\right)\phi_{m'}(\boldsymbol{r} - \boldsymbol{R}')


\\
&\equiv
\sum_{m',\boldsymbol{R}'}c_{m',\boldsymbol{R}'}^{m,\boldsymbol{R}}\phi_{m'}(\boldsymbol{r}-\boldsymbol{R}').
\end{align*}
$$

として展開できます。

この展開係数は、$\boldsymbol{R}' = \boldsymbol{R}$の場合は先ほどと同様に

$$
c_{m',\boldsymbol{R}}^{m,\boldsymbol{R}} \simeq
\varepsilon_m^{\rm a}\delta_{mm'} +\Delta\varepsilon_{m'm} 
$$

で、$\boldsymbol{R}'\neq\boldsymbol{R}$の場合は、同様に「3中心積分」をゼロと置いたりなどして、


$$
c_{m',\boldsymbol{R}'\neq \boldsymbol{R}}^{m,\boldsymbol{R}} \simeq

\int\phi_{m'}^*(\boldsymbol{r} - \boldsymbol{R}')
  V(\boldsymbol{r}-\boldsymbol{R}')
\phi_m(\boldsymbol{r} - \boldsymbol{R})d \boldsymbol{r} 
$$

となります。ここで$\boldsymbol{R}' = \boldsymbol{0}$、つまり「格子点$\boldsymbol{R}$に局在した原子軌道が、原点$\boldsymbol{R}'=\boldsymbol{0}$に飛び移って**行く**確率」を求めてみると、

$$
c_{m',\boldsymbol{0}}^{m,\boldsymbol{R}(\neq \boldsymbol{0})} \simeq

\int\phi_{m'}^*(\boldsymbol{r})
  V(\boldsymbol{r})
\phi_m(\boldsymbol{r} - \boldsymbol{R})d \boldsymbol{r} 
$$

となります。これは先ほど求めた、「原点$\boldsymbol{R} = 0$に局在した原子軌道が、格子点$\boldsymbol{R}'$に飛び移っていく確率$\left(t_{\boldsymbol{R}'}^{m,m'}\right)^*$：


$$
\begin{align*}
-\left(t_{\boldsymbol{R}'}^{m,m'}\right)^* &= 

\left(
\int_V\phi_{m}^*(\boldsymbol{r})
  V(\boldsymbol{r}-\boldsymbol{R}')
\phi_{m'}(\boldsymbol{r}-\boldsymbol{R}')d \boldsymbol{r} 
\right)^*\\

&=

\int_V\phi_{m'}^*(\boldsymbol{r} - \boldsymbol{R}')
  V(\boldsymbol{r})
\phi_m(\boldsymbol{r})d \boldsymbol{r} 
\end{align*}
$$

の複素共役、あるいは複素共役を取る前と同じ形になっています。従って、（複素共役を取る前の）飛び移り積分は「原点に格子点$\boldsymbol{R}$**から**飛び移って**来る**確率」、その複素共役は「原点**から**格子点$\boldsymbol{R}$**に**飛び移って**行く**確率」という意味であったことが分かり、飛び移り積分の複素共役を取る行為は飛び移り方向を逆転させることに対応していたことが分かりました。


## まとめの文章



$$

\int\phi_n^*(\boldsymbol{r} - \boldsymbol{R}_2)V(\boldsymbol{r} - \boldsymbol{R}_1)\phi_m(\boldsymbol{r}-\boldsymbol{R}_1)d\boldsymbol{r}
=
\left(
\int\phi_m^*(\boldsymbol{r} - \boldsymbol{R}_1)V(\boldsymbol{r} - \boldsymbol{R}_1)\phi_n(\boldsymbol{r}-\boldsymbol{R}_2)d\boldsymbol{r}
\right)^*

$$

つまり、格子点$\boldsymbol{R}_{1}$にいる状態$m$の電子が、（微小時間後に）格子点$\boldsymbol{R}_{2} = \boldsymbol{R}_{2} - \boldsymbol{R}$に状態$n$になって飛び移って**来る**確率（に比例（？）する量）を表す

### 飛び移り積分の複素共役の物理的意味

一方、上記飛び移り積分