import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.preprocessing import LabelEncoder
import math
import autograd.numpy as anp
import pymanopt
import pymanopt.manifolds
import pymanopt.solvers

@pymanopt.function.Autograd(manifold=manifold)