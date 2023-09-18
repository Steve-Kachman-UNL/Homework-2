# Homework 2

In the R directory I have included functions for Newton's method and finding a step size that satisfies Wolfe's condition.

### Problems:

1.  On page 53 of the notes. It states that: "Provided $\mathbf{x}^{(t+1)}$ was selected to satisfy the Wolfe's conditions, then $-\mathbf{H}^{(t)}$ being positive definite implies that $-\mathbf{H}^{(t+1)}$ is also positive definite.''

    Show that the above statement is true.

2.  Using the data and model from problem 2.3.

    i.  Derive the log likelihood given in 2.3 a). Note: Censored values use a pmf while uncensored values use a pdf when forming the likelihood function.

    ii. Derive the score function and Hessian for $\alpha$ and $\boldsymbol\beta$.

    iii. Create a function that returns the log likelihood, score function, and Hessian given $\alpha$, ${\boldsymbol\beta}$, $\mathbf{X}$, $\mathbf{y}$, and $\mathbf{w}$.

    iv. Implement a function in R that will compute the MLE given a function that returns the log likelihood, score function, and Hessian of the parameters. The function must also return the standard errors of the parameters.

    v.  Use the function you created to find the MLE and standard errors of $\alpha$ and $\boldsymbol\beta$.

    vi. What do you conclude about the effectiveness of the treatment? How did you arrive at that conclusion?

3.  Implement the Gauss-Newton method for nonlinear regression (section 2.2.3) as an R function.

    -   In addition to estimating the parameters it should also estimate the residual variance.

4.  Using the data from problem 2.6 in the book and the log normal model described in part c.

    i.  Reasonable initial values for $N_0$ and $K$ can be obtained by observing at what $f(t|N_o,K,r)$ reduces to at $t=0$ and $t=\infty$. A reasonable initial value for $r$ can be obtained using equation 2.6.3 and the observed initial growth rate. What are your initial values for $N_0, K, \text{and } r?$ How did you calculate those values?

    ii. Plot $f(t)$ and and the observed data for $t\in [0,200]$ using your initial values. Your plot should be properly labeled.

    iii. Use both Newton's and Gauss-Newton methods to find the MLE of $N_0, K, r, \text{and } \sigma^2.$

    iv. Provide standard errors and estimated correlations for your parameter estimates.
