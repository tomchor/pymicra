Notation
========

Pymicra uses a specific notation to name each one of its columns. This notation
is extremely important, because it is by these labels that Pymicra knows which
variable is in each column, so we advise at the very least a quick look at this chapter to
understand the basics of naming.

The notation can be checked by opening the IPython terminal (because of its
autocomplete function) and entering

\begin{terminal}
\begin{alltt}
import pymicra
dir(pymicra.notation)
\end{alltt}
\end{terminal}

You can change Pymicra's notation at any time by altering the attributes of
\texttt{pymicra.notation}. Some examples are



