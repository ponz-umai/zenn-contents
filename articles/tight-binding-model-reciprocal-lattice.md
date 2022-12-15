---
title: "逆格子ベクトル・逆格子点の物理的な意味"
emoji: "💎"
type: "tech" # tech: 技術記事 / idea: アイデア
topics: ["quantum","quantumcomputing","quantumcomputer","物理","物理学"]
published: true
---
# はじめに

[前章](https://zenn.dev/ponzumai/articles/tight-binding-model-lattice-vec)では結晶の周期を持つ関数のFourier展開について考えました。

結晶が持つ周期は、周期の最小単位を$\boldsymbol{a}_i$として、それらの様々な和

$$
\begin{align*}
\boldsymbol{n} &= n_1\boldsymbol{a}_1 + n_2\boldsymbol{a}_2 + n_3\boldsymbol{a}_3\\
\end{align*}
$$

とする格子ベクトルによって表され、結晶と同じ周期性を持つ関数

$$
V(\boldsymbol{r} + \boldsymbol{n}) = V(\boldsymbol{r})
$$

のFourier級数展開は、逆格子ベクトル$\boldsymbol{k}$を

$$
\begin{align*}
    \boldsymbol{k} 

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

と定義することで、Fourier級数展開を考えることができました。前章では人工的というか、「こうすればうまくいくでしょ」みたいな感じで導入した逆格子ベクトルですが、実際は物理的な意味があるようです。

具体的には、ある結晶の周期性、つまりどのような格子ベクトルで表されるのかを知りたかったとします。その際直接その結晶の格子点を測定するのではなく、実験では逆格子を測定し、そこから逆に格子点を求めることで結晶の格子ベクトルを知ることができます。

# 逆格子の逆格子は直接格子

上述のように、逆格子が得られた場合、「逆格子の逆格子」を求めることで元の結晶の格子点を知ることができます。（結晶の格子点のことを直接格子と呼んだりもします。）

これは直接格子を表す基本併進ベクトルを$\boldsymbol{a}_1,\boldsymbol{a}_2,\boldsymbol{a}_3$に対する
逆格子ベクトル$\boldsymbol{k}$

$$
\begin{align*}
    \boldsymbol{k} 

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

を考え、この$m_1\boldsymbol{b}_1 + m_2\boldsymbol{b}_2 + m_3\boldsymbol{b}_3$を格子ベクトルとみなして、その逆格子を$\boldsymbol{c}_1,\boldsymbol{c}_2,\boldsymbol{c}_3$とすると、

$$
\boldsymbol{c}_1 = \frac{\boldsymbol{b}_2\times \boldsymbol{b}_3}{\boldsymbol{b}_1\cdot(\boldsymbol{b}_2\times\boldsymbol{b}_3)}
$$

などとなります。これは

$$
\boldsymbol{b}_3 = \frac{2\pi\boldsymbol{a}_1\times\boldsymbol{a}_2}{\boldsymbol{a}_3\cdot(\boldsymbol{a}_1\times\boldsymbol{a}_2)}
$$

を代入し、外積の公式

$$
\boldsymbol{A}\times(\boldsymbol{B}\times\boldsymbol{C})
=
\boldsymbol{B}(\boldsymbol{A}\cdot\boldsymbol{C}) - 
\boldsymbol{C}(\boldsymbol{A}\cdot\boldsymbol{B})
$$

を用いると、

$$
\begin{align*}
2\pi\boldsymbol{b}_2\times\boldsymbol{b}_3 
&=
\boldsymbol{b}_2\times \frac{(2\pi)^2\boldsymbol{a}_1\times\boldsymbol{a}_2}{\boldsymbol{a}_3\cdot(\boldsymbol{a}_1\times\boldsymbol{a}_2)}\\

&=
\frac{(2\pi)^2}{\boldsymbol{a}_3\cdot(\boldsymbol{a}_1\times\boldsymbol{a}_2)}
\left\{
    \boldsymbol{a}_1(\boldsymbol{b}_2\cdot \boldsymbol{a}_2) -\boldsymbol{a}_2(\boldsymbol{b}_2\cdot\boldsymbol{a}_1)
\right\}\\

&=
\frac{(2\pi)^3}{\boldsymbol{a}_3\cdot(\boldsymbol{a}_1\times\boldsymbol{a}_2)}
    \boldsymbol{a}_1
\end{align*}
$$

また、上記結果を利用して

$$
\begin{align*}
\boldsymbol{b}_1\cdot(\boldsymbol{b}_2\times\boldsymbol{b}_3)
&=
\frac{(2\pi)^2}{\boldsymbol{a}_3\cdot(\boldsymbol{a}_1\times\boldsymbol{a}_2)}
    \boldsymbol{b}_1\cdot\boldsymbol{a}_1\\
&=\frac{(2\pi)^3}{\boldsymbol{a}_3\cdot(\boldsymbol{a}_1\times\boldsymbol{a}_2)}
\end{align*}
$$

より、結局「逆格子の逆格子」$\boldsymbol{c}_1$は


$$
\boldsymbol{c}_1 = \frac{\boldsymbol{b}_2\times \boldsymbol{b}_3}{\boldsymbol{b}_1\cdot(\boldsymbol{b}_2\times\boldsymbol{b}_3)} = \boldsymbol{a}_1
$$

と、元の直接格子ベクトル$\boldsymbol{a}_1$と一致することが示せます。その他$\boldsymbol{b}_2, \boldsymbol{b}_3$の逆格子も、元の直接格子と一致することが示せます。


というわけである結晶の逆格子を測定することができれば、逆格子の逆格子を求めることで直接格子を知ることができます。
それでは逆格子はどのように測定するのでしょうか。

# 逆格子の測定

そもそも、結晶の直接格子を測定するには、何かわかりませんが超高精度の顕微鏡的なものを利用して原子核や分子・イオンなどを測定し、その周期性を測定する必要があります。しかも3次元結晶であれば表面だけではなく結晶内部まで測定する必要があります。
実験についてはあまり詳しくないですが、まあ、難しそうですよね。

一方、以下のように結晶に色々な波数の光を当てて、入射光と反射光の波数ベクトルを記録することで、対象の結晶の逆格子を知ることができます。

## 結晶によるX線回析のフォン・ラウエによる定式化^[アシュクロフト・マーミン 上I第6章を参照]

以下の黒丸●のように格子点が並んでおり、そこに波数ベクトル$\boldsymbol{k}$の光（X線）を入射します。

この時着目している格子点に対して、あらゆる格子点は格子ベクトルで結ばれています。特に代表して、図では左上の格子点との間の格子ベクトルを$\boldsymbol{n}$とします。

![](/images/tb/r-lattice.png)

入射光は格子点であらゆる方向に散乱されるとして、散乱後に強め合う場合の波数$\boldsymbol{k}'$の条件を考えると、ベクトル$\boldsymbol{n}$で結ばれた格子間で散乱された光同士の行路差が、波長$\lambda = 2\pi/k$の整数倍

$$
n\cos\theta + n\cos\theta' = m\frac{2\pi}{k}
$$

である必要があります。（散乱前後で波数ベクトルの大きさ$|\boldsymbol{k}|, |\boldsymbol{k}'|$は変わらないとします。）

上式は、$kn\cos\theta = -\boldsymbol{k}\cdot\boldsymbol{n},kd\cos\theta' =  \boldsymbol{k}'\cdot\boldsymbol{n}$より、

$$
\boldsymbol{n}\cdot(\boldsymbol{k}' - \boldsymbol{k}) = 2\pi m
$$

となります。上式があらゆる格子ベクトルに対して成り立つという条件は、波数ベクトルの変化

$$
\boldsymbol{K} = \boldsymbol{k}' - \boldsymbol{k}
$$

が逆格子ベクトルと一致するとき、入射光と散乱光が強め合い反射されることを意味しています。

したがって、未知の結晶に対して様々な$\boldsymbol{k}$で光を入射し、入射光に対して反射された光（散乱光同士が強め合った場合の光）の波長を測定し、それらの変化を記録することで、その結晶の直接格子に対する逆格子ベクトルの組を知ることができます。

そうすれば、あとは得られた逆格子ベクトルの逆格子を求めることで結晶の直接格子を知ることができるわけです。

実際の実験ではもう少し色々工夫がされるみたいですが、今のところはいったんこれくらいにとどめて置いて、逆格子への親近感を高められたところで終わりにします。

# おわりに

今回の内容は、モデルの理解という意味では直接必要な内容ではなかったようにも思うのですが、逆格子（ベクトル）について分かったような分らんような感覚をずっと持っていたので折角ならちゃんと理解しておこうと、書いてみました。

とはいえ現状でも誤魔化していたり、深くは入らなかったりしているのですが、ひとまずはこんなところで次に進んでいこうと思います。

次章では周期的なポテンシャルがある場合の電子の性質についてまとめる予定です。