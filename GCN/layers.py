from mxnet.gluon import HybridBlock
from mxnet.gluon.nn import Activation
import mxnet.ndarray as nd
import mxnet as mx

class SpectralRule(HybridBlock):
    def __init__(self, A, in_units, out_units, activation='relu', **kwargs):
        super().__init__(**kwargs)
# no need to introduce self-loops by adding A to identity matrix. every cell's interaction value with itself is 
# already included in A
# no need to normalize by inverse degree matrix, as values in adjacency matrix are already normalized 
        A_hat = A.copy()
#         D = np.array(np.sum(A_hat, axis=0))
    #         D = nd.sum(A_hat, axis=0)
#         D_inv = D**-0.5
#         D_inv = np.matrix(np.diag(D_inv))
#         D_inv = nd.diag(D_inv)
        
#         D_hat = np.array(np.sum(A, axis=0))
#         D_hat = np.matrix(np.diag(D_hat))

#         A_hat = D_inv * A_hat * D_inv
        
        self.in_units, self.out_units = in_units, out_units
        
        with self.name_scope():
            self.A_hat = self.params.get_constant('A_hat', A_hat)
            self.W = self.params.get(
                'W', shape=(self.in_units, self.out_units)
            )
            if activation == 'identity':
                self.activation = lambda X: X
            else:
                self.activation = Activation(activation)

    def hybrid_forward(self, F, X, A_hat, W):
        aggregate = F.dot(A_hat, X)
        propagate = self.activation(
            F.dot(aggregate, W))
        return propagate