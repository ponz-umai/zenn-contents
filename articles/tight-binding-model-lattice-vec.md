---
title: "格子ベクトルと逆格子ベクトル"
emoji: "🐥"
type: "tech" # tech: 技術記事 / idea: アイデア
topics: ["quantum","quantumcomputing","quantumcomputer","物理","物理学"]
published: false
---
# はじめに
前回記事で、tight-binding近似の基本となる原子核に強く束縛された多電子状態の波動関数について理解しました。これで早速、固体の中の電子状態に入っていけます。
そこで早速固体中の電子が取る状態について進んでいきます。

手始めに今回は数学的な準備として、固体、特に周期的な結晶が持つ周期性の取り扱いについて記載していきます。具体的には、結晶中の物理量やハミルトニアンは、結晶と同じ周期性を持っているはずです。そのような関数のFourier変換について整理します。

# 固体の中の周期性





# Bravais格子ベクトルと逆格子ベクトル
本稿で最終的なゴールとしている固体は、ほとんどの場合原子や分子、イオンが周期的に配列した状態を想定しています。
それらの座標を**格子点**等と呼び、また座標に対応する位置ベクトルを **(Bravais)格子ベクトル**とよびます。

各格子ベクトルは一定の方向に周期性を持つので、最小単位となるベクトル$\boldsymbol{a}_i$を用いて

$$
\boldsymbol{n} = n_1\boldsymbol{a}_1 + n_2\boldsymbol{a}_2 + n_3\boldsymbol{a}_3
$$

のように書けます。このような各方向への周期性に対応するベクトル$\boldsymbol{a}_i$を基本並進ベクトルと呼び、また、全ての$n_i$の集合をBraveis格子と呼び、固体の周期性（多くの場合対称性と言ったりもします）を表す指標になります。

まあ、そういう細かい名称はしっかりした教科書に譲ることにして、本題に行きましょう。


## 物理量の周期性

固体が上記のような周期性を持っていることから、固体中の物理量、例えば固体中の電子が感じるポテンシャル$V(\boldsymbol{r})$は、固体の周期性と同じく

$$
V(\boldsymbol{r} + \boldsymbol{n}) = V(\boldsymbol{r})
$$

のような周期関数となっているはずです。このような関数のFourier級数展開が今回の主題です。




## 1次元関数のFourier級数展開

初めにもっとも簡単な例として

$$
f(x + a) = f(x)
$$

のような1次元の周期のみを考えると、これはよくあるFourier級数展開を用いて

$$
f(x) = \sum_k C_k e^{ikx},\\
k = \frac{2\pi}{a}n, n = 0,1,2,\cdots,\\
C_k = \frac{1}{a}\int_0^a f(x)e^{-ikx}dx 
$$

と表せます^[なお、Fourier級数展開が元の関数と一致するこの証明については[こちらのページ](http://www.maroon.dti.ne.jp/koten-kairo/works/fft/fft_start.html)や、[こちらのページ](http://www.core.kochi-tech.ac.jp/m_inoue/work/pdf/sekiguti/colleage/3.pdf)に詳しい解説があります。]


## 結晶格子の周期を持つ関数のFourier 級数展開

### いくつかの基本並進ベクトル

天下り的ですが、代表的な結晶構造に対して以下のような基本併進ベクトルが対応しています。なお以下で$\boldsymbol{e}_i$は$i = x,y,z$とし、それぞれ直交座標の$x,y,z$軸に対する単位ベクトルを表します。

#### 単純立方格子

$$
\boldsymbol{a}_1 = a\boldsymbol{e}_x, \\
\boldsymbol{a}_2 = a\boldsymbol{e}_y,\\
\boldsymbol{a}_3 = a\boldsymbol{e}_z.
$$

これは最も簡単な場合ですね。本稿の範囲では多分、結局ほとんどを単純立方格子の例で済ませることになるとは思います。

#### 面心立方格子

$$
\boldsymbol{a}_1 = \frac{a}{2}\boldsymbol{e}_x + \frac{a}{2}\boldsymbol{e}_y, \\ 
{}\\
\boldsymbol{a}_2 =  \frac{a}{2}\boldsymbol{e}_y + \frac{a}{2}\boldsymbol{e}_z,\\
{}\\
\boldsymbol{a}_3 =  \frac{a}{2}\boldsymbol{e}_x + \frac{a}{2}\boldsymbol{e}_z.
$$

#### 体心立方格子

$$
\boldsymbol{a}_1 = \frac{a}{2}\boldsymbol{e}_x + \frac{a}{2}\boldsymbol{e}_y
+ \frac{-a}{2}\boldsymbol{e}_z\\ 
{}\\
\boldsymbol{a}_2 =  \frac{-a}{2}\boldsymbol{e}_x + \frac{a}{2}\boldsymbol{e}_y
+ \frac{a}{2}\boldsymbol{e}_z,\\
{}\\
\boldsymbol{a}_3 =  \frac{a}{2}\boldsymbol{e}_x + \frac{-a}{2}\boldsymbol{e}_y
+ \frac{a}{2}\boldsymbol{e}_z.
$$

### 単純立方格子の場合のFourier級数展開

初めに最も簡単な例として、単純立方格子の基本併進ベクトルで表される周期をもつ関数$V(\boldsymbol{r})$,つまり格子ベクトルを

$$
\begin{align*}
\boldsymbol{n} &= n_1\boldsymbol{a}_1 + n_2\boldsymbol{a}_2 + n_3\boldsymbol{a}_3\\
&=n_1a\boldsymbol{e}_x + n_2a\boldsymbol{e}_y + n_3a\boldsymbol{e}_z
\end{align*}
$$

としたとき任意の$\boldsymbol{n}$について

$$
V(\boldsymbol{r} + \boldsymbol{n}) = V(\boldsymbol{r})
$$

となる関数ののFourier展開を考えます。

$\boldsymbol{r}$を正規直交基底$\boldsymbol{e}_i$で成分表示すると、

$$
V(x + n_1a,y,z) =V(x,y+n_2a,z) = V(x,y,z+n_3a) = V(x,y,z) 
$$

と、それぞれの座標成分に対して周期性を持ちます。そこで順次$x$座標からFourier展開していきます。

以下で$k_i= am_i/2\pi, m_i = 0,1,2\cdots (i = x,y,z)$として、

$$
\begin{align*}
V(x,y,z) &= \sum_{k_x} C_{k_x}(y,z)e^{ik_xx}\\
C_{k_x}(y,z) &= \frac{1}{a}\int_0^aV(x,y,z)dx
\end{align*}
$$

さらに$C_{k_x}(y+n_2a,z) = C_{k_x}(y,z)$なので、


$$
\begin{align*}
C_{k_x}(y,z) &=\sum_{k_y}C_{k_x,k_y}(z)e^{ik_yy},\\
C_{k_x,k_y}(z) &= \frac{1}{a}\int_0^a C_{k_x}(y,z)e^{-ik_yy}dy\\

&= \frac{1}{a^2}\int_0^a \int_0^a V(x,y,z)e^{-i(k_xx + k_yy)}dxdy

\end{align*}
$$

さらにさらに（略）

$$
\begin{align*}
C_{k_x,k_y}(z) &=\sum_{k_z}C_{k_x,k_y,k_z}e^{ik_zz},\\
C_{k_x,k_y,k_z} &= \frac{1}{a}\int_0^a C_{k_x,k_y}(z)e^{-ik_zz}dz\\

&= \frac{1}{a^3}\int_0^a\int_0^a\int_0^a V(x,y,z)e^{-i(k_xx + k_yy + k_zz)}dxdydz
\end{align*}
$$

全部合わせて

$$
V(x,y,z) = \sum_{k_x,k_y,k_z}C_{k_x,k_y,k_z}e^{i(k_xx + k_yy + k_zz)}\\
C_{k_x,k_y,k_z}= \frac{1}{a^3}\int_0^a\int_0^a\int_0^a V(x,y,z)e^{-i(k_xx + k_yy + k_zz)}dxdydz
$$

となります。

ただしこのままだと位置ベクトル$\boldsymbol{r}$に対して、基底を選んで成分表示する必要がある形になっています。そこで

$$

\boldsymbol{k} = k_x\boldsymbol{e}_x + k_y\boldsymbol{e}_y + k_z\boldsymbol{e}_z

$$

というベクトルを定義すると、

$$
\begin{align*}
\boldsymbol{k}\cdot\boldsymbol{r} &= k_x\boldsymbol{e}_x\cdot\boldsymbol{r} + 
k_y\boldsymbol{e}_y\cdot\boldsymbol{r} + 
k_z\boldsymbol{e}_z\cdot\boldsymbol{r} \\
&= k_xx + k_yy + k_zz
\end{align*}
$$

と基底を気にせずにベクトルの内積の形で書くことができます。以上より単純立方格子の周期をもつ関数のFourier展開は

$$
V(\boldsymbol{r}) = \sum_{\boldsymbol{k}}C_{\boldsymbol{k}}e^{i\boldsymbol{k}\cdot\boldsymbol{r}}\\
C_{\boldsymbol{k}}= \frac{1}{v_c}\int_{v_c} V(\boldsymbol{r})e^{-i\boldsymbol{k}\cdot\boldsymbol{r}}d\boldsymbol{r}
$$

と、シンプルな形で書けました。ただし、$\sum_{\boldsymbol{k}}$は$k_x,k_y,k_z$それぞれの総和、$v_c$は繰り返しの最小単位の体積を表し、$v_c = a^3, \int_{v_c}d\boldsymbol{r}=\int_0^adx\int_0^ady\int_0^adz$です。^[物性でよく使われる記号なのでここでも採用しています。$v_c$の$c$は、多分最小単位の英語表記"unit cell"のcだと思います]

### 非直交な周期性を持つ関数のFourier展開(2次元）

ここまでは特に引っかかることもなく進めましたが、面心立方格子や体心立方格子のような周期性を扱うには、直交しない方向への周期性を考える必要があります。

まずは簡単に、2次元の場合を考えてみます。

関数$V(\boldsymbol{r})$が図のような、

![](/images/tb/fourier1.png)

基本併進ベクトル
$\boldsymbol{a}_1 , \boldsymbol{a}_2$に対して周期性を持つ場合のFourier展開について考えてみます。

ここで$|\boldsymbol{a}_1| =a_1, |\boldsymbol{a}_2|=a_2$として関数の周期と平行な単位ベクトル

$$
\boldsymbol{e}_1 = \frac{\boldsymbol{a}_1}{a_1}, \boldsymbol{e}_2 = \frac{\boldsymbol{a}_2}{a_2}
$$

を基底として成分表示をすることにします。

$$
\begin{align*}
\boldsymbol{r} &= r_1\boldsymbol{e}_1 + r_2\boldsymbol{e}_2\\

&= (r_1,r_2)
\end{align*}
$$

このような斜交基底を定義しても任意の位置$\boldsymbol{r}$を表すことができますが、$r_i$は、直交基底の場合と異なり$\boldsymbol{e}_1\cdot\boldsymbol{r}$ではないことに注意してください。（じゃあどうすればいいのかというと、後でわかります）

![](/images/tb/fourier2.png)


このように成分表示すると、周期性は

$$
V(r_1+n_1a_1,r_2) = V(r_1,r_2 + n_2a_2) = V(r_1,r_2)
$$

と書けます。そこで先ほどと同じように各成分に対してFourier展開すると、

$$
V(r_1,r_2) = \sum_{k_1,k_2}C_{k_1,k_2}e^{i(k_1r_1 + k_2r_2)}\\
C_{k_1,k_2}= \frac{1}{a_1a_2}\int_0^{a_1}\int_0^{a_2} V(r_1,r_2)e^{-i(k_1r_1 + k_2r_2)}dr_1dr_2,\\
k_i = \frac{2\pi}{a_i}m_i, m_i = 0,1,2,\cdots
$$

と書けます。ところで今回も基底ベクトルに依存した成分表示になっており、特に斜交座標の成分表示なので使いにくい形になっています。

そこで先ほどと同じように、基底座標を意識しなくても自動的に上記のような形になるようなベクトル$\boldsymbol{k}$を定義したいのですが、どうすれば良いでしょうか。

素朴に

$$
\widetilde{\boldsymbol{k}}= k_1\boldsymbol{e}_1+k_2\boldsymbol{e}_2
$$

としてみても、今回は先述のように$\boldsymbol{e}_1\cdot\boldsymbol{r}=r_1$とできず、

$$
\widetilde{\boldsymbol{k}}\cdot\boldsymbol{r} \neq k_1r_1 + k_2r_2
$$

となってしまいます。そこで上手いこと、

$$
\boldsymbol{e}_i^* \cdot\boldsymbol{r} = r_i
$$


となってくれるようなベクトル$\boldsymbol{e}_i^*$を見つければ、これを用いて

$$
\boldsymbol{k} = k_1\boldsymbol{e}_1^* + k_2\boldsymbol{e}_2^*
$$

と置けば、

$$
\boldsymbol{k} \cdot \boldsymbol{r}

=
k_1\boldsymbol{e}_1^*\cdot\boldsymbol{r}
+
k_2\boldsymbol{e}_2^*\cdot\boldsymbol{r}
=
k_1r_1 + k_2r_2
$$

と直交座標の場合と同じ形に書けそうです。



そこで、$\boldsymbol{e}_i^* \cdot\boldsymbol{r} = r_i$を満たす$\boldsymbol{e}_i^*$を考えるために改めて$\boldsymbol{r}$を$\boldsymbol{e}_1,\boldsymbol{e}_2$で展開してみると、

$$
\boldsymbol{r} = r_1\boldsymbol{e}_1 + r_2\boldsymbol{e}_2
$$

なので、$\boldsymbol{e}_1$との内積が$1$で、$\boldsymbol{e}_2$に直交するようなベクトル$\boldsymbol{e}^*_1$と、逆に
$\boldsymbol{e}_2$との内積が$1$で、$\boldsymbol{e}_1$に直交するベクトル$\boldsymbol{e}^*_2$を考えればよさそうです。

まず$\boldsymbol{e}_2$と直交するという条件を考えてみると、「2つのベクトルと直交する方向を持つベクトル」として定義される外積が使えそうです。
今回は$\boldsymbol{e}_2$と直交してさえいればいいので、もう一つのベクトルの向きは（$\boldsymbol{e}_1$と平行じゃなければ）どこでもいいのですが、$\boldsymbol{r}$と同じ平面上にある方が良さそうなので$\boldsymbol{e}_1, \boldsymbol{e}_2$と直交し大きさが$1$のベクトル$\boldsymbol{n}$を導入し、$A$を後で調整する定数として

$$
\boldsymbol{e}_1^* = A\boldsymbol{e}_2\times\boldsymbol{n}

$$

$\boldsymbol{e}_2^*$についても同様に

$$
\boldsymbol{e}_2^* = B\boldsymbol{e}_1\times\boldsymbol{n}
$$

と置けば良さそうです。

続いて$\boldsymbol{e}_1^*\cdot\boldsymbol{e}_1=\boldsymbol{e}_2^*\cdot\boldsymbol{e}_2=1$の条件を考えると、

$$
\boldsymbol{e}_1^*\cdot\boldsymbol{e}_1
=
A(\boldsymbol{e}_2\times\boldsymbol{n})
\cdot \boldsymbol{e}_1 = 1\\
\boldsymbol{e}_2^*\cdot\boldsymbol{e}_2
=
B(\boldsymbol{e}_1\times\boldsymbol{n})
\cdot \boldsymbol{e}_2 = 1
$$

を満たせばよいので、$(\boldsymbol{e}_i\times\boldsymbol{n})\cdot \boldsymbol{e}_j$はスカラーなのでその逆数を選べば良くて、最終的に

$$
\boldsymbol{e}_1^* =  \frac{\boldsymbol{e}_2\times\boldsymbol{n}}{\boldsymbol{e}_1\cdot(\boldsymbol{e}_2\times\boldsymbol{n})}, 

\boldsymbol{e}_2^* =  \frac{\boldsymbol{e}_1\times\boldsymbol{n}}{\boldsymbol{e}_2\cdot(\boldsymbol{e}_1\times\boldsymbol{n})},
$$

と定義すれば、任意のベクトル$\boldsymbol{r}$に対して

$$
\boldsymbol{e}_1^*\cdot\boldsymbol{r}=r_1,\\

\boldsymbol{e}_2^*\cdot\boldsymbol{r}=r_2
$$

が満たされることになります。

したがって、

$$
\begin{align*}
\boldsymbol{k} = 
k_1\boldsymbol{e}_1^*
+
k_2\boldsymbol{e}_2^*
=

k_1\frac{\boldsymbol{e}_2\times\boldsymbol{n}}{\boldsymbol{e}_1\cdot(\boldsymbol{e}_2\times\boldsymbol{n})}

+
k_2\frac{\boldsymbol{e}_1\times\boldsymbol{n}}{\boldsymbol{e}_2\cdot(\boldsymbol{e}_1\times\boldsymbol{n})}
\end{align*}
$$

と定義すれば^[このように定義したベクトルを「双対ベクトル」等と呼ぶようです。]、
$$
\boldsymbol{k}\cdot\boldsymbol{r} = k_1r_1 + k_2r_2
$$
を満たし、

周期$\boldsymbol{a}_1,\boldsymbol{a}_2$を持つ2次元関数

$$
V(\boldsymbol{r} + \boldsymbol{a}_1) 

= V(\boldsymbol{r} + \boldsymbol{a}_2)

= V(\boldsymbol{r})
$$

に対して、Fourier展開を

$$
\begin{align*}
V(\boldsymbol{r}) &= \sum_{k_1,k_2}C_{k_1,k_2}e^{i(k_1r_1 + k_2r_2 )}
=
\sum_{\boldsymbol{k}}C_{\boldsymbol{k}}e^{i\boldsymbol{k}\cdot\boldsymbol{r}}\\

C_{k_1,k_2}&= \frac{1}{a_1a_2}\int_0^{a_1}\int_0^{a_2} V(r_1,r_2)e^{-i(k_1r_1 + k_2r_2 )}dr_1dr_2\\

&=
\frac{1}{s_c}\int_{s_c} V(\boldsymbol{r})e^{-i\boldsymbol{k}\cdot\boldsymbol{r}}d\boldsymbol{r}
\end{align*}
$$

と書けることがわかりました。
なお、ここでは斜交座標なので$\boldsymbol{a}_1$と$\boldsymbol{a}_2$のなす角を$\theta$として、繰り返しの最小単位の面積を$s_c = a_1a_2\sin\theta$とし、
積分
$\int_{s_c}d\boldsymbol{r}=\int_0^{a_1}\int_0^{a_2}\sin\theta dr_1dr_2$です。


### 非直交な周期性を持つ関数のFourier展開(3次元）

2次元の場合と同様に、3次元の非直交な周期
$\boldsymbol{a}_1, \boldsymbol{a}_2, \boldsymbol{a}_3$
を持つ関数

$$
V(\boldsymbol{r}+\boldsymbol{a}_1) 
= 
V(\boldsymbol{r}+\boldsymbol{a}_2) 
= 
V(\boldsymbol{r}+\boldsymbol{a}_3) 
= 
V(\boldsymbol{r})
$$

も、$\boldsymbol{a}_i$方向の単位ベクトルを$\boldsymbol{e}_i$として、それらを基底に

$$
\boldsymbol{r} = \sum_ir_i\boldsymbol{e}_i
$$

と置くと、

「$\boldsymbol{e}_2$と$\boldsymbol{e}_3$に直交して$\boldsymbol{e}_1$との内積が$1$のベクトル」$\boldsymbol{e}_1^*$、

「$\boldsymbol{e}_1$と$\boldsymbol{e}_3$に直交して$\boldsymbol{e}_2$との内積が$1$のベクトル」$\boldsymbol{e}_2^*$、

「$\boldsymbol{e}_1$と$\boldsymbol{e}_2$に直交して$\boldsymbol{e}_3$との内積が$1$のベクトル」$\boldsymbol{e}_3^*$

を用いて、

$$
\boldsymbol{k} = k_1\boldsymbol{e}_1^*
+

k_2\boldsymbol{e}_2^*
+
k_3\boldsymbol{e}_3^*
$$

として、Fourier展開

$$
V(\boldsymbol{r}) = \sum_{\boldsymbol{k}}C_{\boldsymbol{k}}e^{i\boldsymbol{k}\cdot\boldsymbol{r}}\\
C_{\boldsymbol{k}}= \frac{1}{v_c}\int_{v_c} V(\boldsymbol{r})e^{-i\boldsymbol{k}\cdot\boldsymbol{r}}d\boldsymbol{r}
$$

と表すことができます。ただし、$\boldsymbol{e}_i^*$はそれぞれ、

$$
\boldsymbol{e}_1^*
=
\frac{\boldsymbol{e}_2\times\boldsymbol{e}_3}{\boldsymbol{e}_1\cdot(\boldsymbol{e}_2\times\boldsymbol{e}_3)}

\boldsymbol{e}_2^*
=
\frac{\boldsymbol{e}_3\times\boldsymbol{e}_1}{\boldsymbol{e}_2\cdot(\boldsymbol{e}_3\times\boldsymbol{e}_1)}

\boldsymbol{e}_3^*
=
\frac{\boldsymbol{e}_1\times\boldsymbol{e}_2}{\boldsymbol{e}_3\cdot(\boldsymbol{e}_1\times\boldsymbol{e}_2)}
$$

で、$v_c$は$\boldsymbol{a}_1,\boldsymbol{a}_2,\boldsymbol{a}_3$で囲まれる平行六面体の体積$\boldsymbol{a}_1\cdot(\boldsymbol{a}_2\times\boldsymbol{a}_3)$で、$1/v_c\int_{v_c}d\boldsymbol{r} = 1/(a_1a_2a_3)\int_0^{a_1}dr_1\int_0^{a_2}dr_2\int_0^{a_3}dr_3$です。

ここで、序盤で考えたように直交した周期性を持つ場合は、

$$
\boldsymbol{e}_2\times\boldsymbol{e}_3 = \boldsymbol{e}_1, \boldsymbol{e}_3\times\boldsymbol{e}_1 = \boldsymbol{e}_2,\boldsymbol{e}_1\times\boldsymbol{e}_2 = \boldsymbol{e}_3
$$

となることから$\boldsymbol{e}_i^* = \boldsymbol{e}_i$となります。

## 逆格子ベクトル

以上で、任意の周期を持つ関数に対してFourier展開を定義することができました。

ここで、

$$
k_i = \frac{2\pi}{a_i}m_i, m_i = 0,1,2,\cdots
$$

と

$$
\frac{\boldsymbol{e}_j\times\boldsymbol{e}_l}{a_i\boldsymbol{e}_i\cdot(\boldsymbol{e}_j\times\boldsymbol{e}_l)}

=
\frac{\boldsymbol{a}_j\times\boldsymbol{a}_l}{\boldsymbol{a}_i\cdot(\boldsymbol{a}_j\times\boldsymbol{a}_l)}
$$

から、$\boldsymbol{k}$を

$$
\begin{align*}
    \boldsymbol{k} 
    &= k_1\frac{\boldsymbol{e}_2\times\boldsymbol{e}_3}{\boldsymbol{e}_1\cdot(\boldsymbol{e}_2\times\boldsymbol{e}_3)}
+

k_2\frac{\boldsymbol{e}_3\times\boldsymbol{e}_1}{\boldsymbol{e}_2\cdot(\boldsymbol{e}_3\times\boldsymbol{e}_1)}
+
k_3\frac{\boldsymbol{e}_1\times\boldsymbol{e}_2}{\boldsymbol{e}_3\cdot(\boldsymbol{e}_1\times\boldsymbol{e}_2)}\\

&=
     m_1\frac{2\pi\boldsymbol{e}_2\times\boldsymbol{e}_3}{a_1\boldsymbol{e}_1\cdot(\boldsymbol{e}_2\times\boldsymbol{e}_3)}
+

 m_2\frac{2\pi\boldsymbol{e}_3\times\boldsymbol{e}_1}{a_2\boldsymbol{e}_2\cdot(\boldsymbol{e}_3\times\boldsymbol{e}_1)}
+
 m_3\frac{2\pi\boldsymbol{e}_1\times\boldsymbol{e}_2}{a_3\boldsymbol{e}_3\cdot(\boldsymbol{e}_1\times\boldsymbol{e}_2)}\\

&=
 m_1\frac{2\pi\boldsymbol{a}_2\times\boldsymbol{a}_3}{\boldsymbol{a}_1\cdot(\boldsymbol{a}_2\times\boldsymbol{a}_3)}
+

 m_2\frac{2\pi\boldsymbol{a}_3\times\boldsymbol{a}_1}{\boldsymbol{a}_2\cdot(\boldsymbol{a}_3\times\boldsymbol{a}_1)}
+
 m_3\frac{2\pi\boldsymbol{a}_1\times\boldsymbol{a}_2}{\boldsymbol{a}_3\cdot(\boldsymbol{a}_1\times\boldsymbol{a}_2)}\\

 &\equiv
 m_1\boldsymbol{b}_1 + m_2\boldsymbol{b}_2 + m_3\boldsymbol{b}_3
\end{align*}
$$

とも書けます。これならいちいち「$\boldsymbol{a}_i$方向の単位ベクトル」等と考える必要もなく、与えられた周期に対応するベクトル$\boldsymbol{a}_i$を使って直接波数ベクトル$\boldsymbol{k}$を定義できます。

物性物理（固体物理）の分野では、結晶の周期を表すベクトル

$$
\boldsymbol{n} = n_1\boldsymbol{a}_1
+
n_2\boldsymbol{a}_2

+
n_3\boldsymbol{a}_3

$$

を格子ベクトルと呼ぶことに対応して、このようなベクトル

$$
\boldsymbol{k} = m_1\boldsymbol{b}_1 + m_2\boldsymbol{b}_2 + m_3\boldsymbol{b}_3
$$

を**逆格子ベクトル**と呼び、$m_i$を指定して得られる各$\boldsymbol{k}$が示す座標を**逆格子点**と呼びます。

また、逆格子ベクトルであることを強調する意味で小文字の$\boldsymbol{k}$ではなく、大文字の$\boldsymbol{K}$や$\boldsymbol{G}$と表すことが多いです。

またこの時任意の格子ベクトル$\boldsymbol{n}$に対して、

$$
\boldsymbol{b}_i\cdot\boldsymbol{a}_j = 2\pi\delta_{ij} 
$$

となることから、

$$
\boldsymbol{K}\cdot\boldsymbol{r} = 2\pi l, l = 0,1,2,\cdots
$$

となります。

## 例題　面心立方格子の結晶内のポテンシャル

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

# まとめ

今回は普通の教科書ではさらっと出てくる逆格子ベクトルについて、（必要以上に？）詳しく考えてみました。

改めてまとめると、



これで