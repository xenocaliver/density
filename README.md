# density evolution demo program

## Preparation
Install following tools.
- gcc, g++
- make
- cmake
- boost
- [Qt6](https://www.qt.io/product/qt6)
- [qwt](https://qwt.sourceforge.io/index.html)
- [fftw3](http://fftw.org)

gcc, g++, make, cmake, boost, Qt6, and qwt are installed by package management system. [fftw3](http://fftw.org) is downloaded and compiled during this program's compilation process.

## Downdloading the demo program
Clone demo program's repository.

```sh
git clone https://github.com/xenocaliver/density.git ./
```
## Build the demo program
Before building the demo program, you must modify CMakeLists.txt file for your system. This CMakeLists.txt is written for mac OS and homebrew.

```sh
cd density
mkdir build
cd build
cmake ..
make
```

## Run the demo program
Run demo program as follows.

```sh
./density <degree distribution polynomial json file> <sigma>
```

## Implimentation details
I implemented [Chung](https://ieeexplore.ieee.org/document/905935)'s denisty evolution method. And this program is fork of Professor [Wadayama](https://wadayama.github.io)'s code.
I added GUI plotting pdf(probability density function) convergence process and adopted to irregular LDPC codes. 

Let $v$ be a log-likelihood ratio (LLR) message from a degree-$d_{v}$ variable node to a check node. Due to BP(Belief Propagation) decoding method, $v$ equals to the sum of all incoming LLRs; i.e.

$$
v=\sum_{i=0}^{d_{v}-1}u_{i}.
$$

And check node update rule is shown as follows:

$$
u = 2\tanh^{-1}\left[\prod_{j=1}^{d_{c}-1}\tanh\frac{v_{j}}{2}\right]
$$

where $v_{j}, j=1,\ldots,d_{c}-1$.

Let $\mathit{\Delta}$ be a quantization width of random variables. And we define quantize operator $\mathcal{Q}$ as follows:

$$
\mathcal{Q}(w)=
\begin{cases}
\left\lfloor\frac{w}{\mathit{\Delta}}+\frac{1}{2}\right\rfloor\mathit{\Delta}& \left(w\geqq\frac{\mathit{\Delta}}{2}\right)\\
\left\lceil\frac{w}{\mathit{\Delta}}-\frac{1}{2}\right\rceil\mathit{\Delta}& \left(w\leqq-\frac{\mathit{\Delta}}{2}\right).
\end{cases}
$$

And let $\bar{w}$ be a quantized message of $w$ i.e. $\bar{w}=\mathcal{Q}(w)$. We denote the probability mass function(pmf) of a quantized message $\bar{w}$ by $p_{w}[k]=\mathrm{Pr}(\bar{w}=k\mathit{\Delta})$ for $k\in\mathbb{Z}$. Then, let $p_{v}$ be a pmf corresponding message $v$ and $p_{v}$ is given by

$$
p_{v}=\bigotimes_{i=0}^{d_{v}-1}p_{u_{i}}
$$

where $\bigotimes$ represents convolution calculation. In general, computing complexity for convolution is very large. However, we used (complex) Fourier transform $\mathscr{F}$. Let $\mathscr{F}[w]$ be a (complex) Fourier transformed function. Therefore, we apply Fourier transform both side of above fomula and we have

$$
\mathscr{F}[p_{v}]=\mathscr{F}\left[\bigotimes_{i=0}^{d_{v}-1}p_{u_{i}}\right]=\prod_{i=0}^{d_{v}-1}\mathscr{F}[p_{u_{i}}].
$$

So, we apply (complex) Fourier transform once, and calculate simple product of transformed functions and revert the result via inverse Fourier transform $\mathscr{F}^{-1}$. Thus, we can calculate variable node update of pdf.

However, check node update calculation is more complicated. Before we formulate update rule for check node update, we define

$$
\mathcal{R}(a,b)=\mathcal{Q}\left(2\tanh^{-1}\left(\tanh\left(\frac{a}{2}\right)\tanh\left(\frac{b}{2}\right)\right)\right).
$$

Then, we can convert check node update rule into 

$$
\bar{u}=\mathcal{R}\left(\bar{v}_{1},\mathcal{R}(\bar{v}_{2},\mathcal{R}(\bar{v}_{3},\ldots,\mathcal{R}(\bar{v}_{d_{c}-2},\bar{v}_{d_{c}-1}))\right).
$$

In order to execute above calculation, let $c=\mathcal{R}(a,b)$ and $\mathcal{R}$ calculation in quantized form:

$$
p_{c}[k] = \sum_{(i,j):k\mathit{\Delta}=\mathcal{R}(i\mathit{\Delta},j\mathit{\Delta})}p_{a}[i]p_{b}[j].
$$

We construct lookup table for $(i,j)$ satisfies $k\mathit{\Delta}=\mathcal{R}(i\mathit{\Delta},j\mathit{\Delta})$ before evolution calculation.

Above all, we execute quantized density evolution of pdfs variable node update and check node update back and forth.

## Irregular LDPC specifications

An irregular LDPC code's profile is specified by degree distribution polynomials. Variable node degree distribution polynomial in edge perspective is defined by

$$
\lambda(x)=\sum_{i}\lambda_{i}x^{i-1}
$$

and check node degree distribution polynomial is defined by

$$
\rho(x)=\sum_{i}\rho_{i}x^{i-1}.
$$

Where $0\leqq\lambda_{i},\rho_{i} < 1$ and 

$$
\begin{aligned}
\lambda(1) &=\sum_{i}\lambda_{i} = 1\\
\rho(1)&=\sum_{i}\rho_{i} = 1.
\end{aligned}
$$

I expressed degree distribution polynomial by a json file. For example, $(3,6)-$ regular LDPC code's degree distribution polynomials are shown as follows:

$$
\begin{aligned}
\lambda(x) &= x^{2}\\
\rho(x)&=x^{5}.
\end{aligned}
$$

And an example of optimized irregular LDPC code's degree distribution polynomials are

$$
\begin{aligned}
\lambda(x)&=0.18411x+0.50891x^{2}+0.09001x^{14}+0.21697x^{15}\\
\rho(x)&=0.5x^{3}+0.5x^{4}.
\end{aligned}
$$

The content of json file corresponding to the first example is shonw below.

```
{
    "degree_distribution":
    {
        "lambda":
            [{"degree": 2, "weight": 1.0}],
        "rho":
            [{"degree": 5, "weight": 1.0}]
    }
}
```

The content of json file corresponding to the second example is shown below.

```
{
    "degree_distribution":
    {
        "lambda":
        [
                {"degree": 1, "weight": 0.18411},
                {"degree": 2, "weight": 0.50891},
                {"degree": 14, "weight": 0.09001},
                {"degree": 15, "weight": 0.21697}
        ],
        "rho":
            [
                {"degree": 3, "weight": 0.5},
                {"degree": 4, "weight": 0.5}
            ]
    }
}
```

This program read these json file and reflect density evolution update rules.