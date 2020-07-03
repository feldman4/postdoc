from rtRosetta.v12_simple import *
from ..sequence import print_alignment

MODEL_PATH = 'wfc/models/model_{token}.npy'

def argmax_grad_bypass(y):
    # returns argmax but gradient comes from softmax
    # output is sum of softmax (w/ gradient) and argmax-softmax
    #   (no gradient)
    y = tf.nn.softmax(y, -1)
    y_hard = tf.one_hot(tf.argmax(y, -1), 20)   # argmax
    return tf.stop_gradient(y_hard-y)+y       # gradient bypass


def setup_model(num_models, loss_fcn, inputs, tokens=[], model_path=MODEL_PATH, 
    clear=True):
    if clear:
        K.clear_session()

    # inputs
    I = Input(shape=(None, 20), dtype=tf.float32)
    
    # forward pass uses the argmax of PSSM representing search state
    I_hard = argmax_grad_bypass(I)
    
    # inserts zeros, batch,length,20 to batch,length,21
    I_hard_gap = tf.pad(I_hard, [[0, 0], [0, 0], [0, 1]])  # add gap

    # load each model
    avg_feat, avg_bb = [], []

    if not tokens:
        tokens = ["xaa", "xab", "xac", "xad", "xae"][:num_models]
    for token in tokens:
        print(f"loading model: {token}")
        weights = load_weights(model_path.format(token=token), mode="TrRosetta")
        feat, bb = RESNET(weights=weights)(I_hard_gap)
        avg_feat.append(feat)
        avg_bb.append(bb)

    feat = tf.reduce_mean(avg_feat, 0)
    bb = K.sum(tf.reduce_mean(avg_bb, 0)[..., 1:], -1)

    loss = loss_fcn(feat)
    
    # compute gradients given sum of loss
    grad = Lambda(lambda x: tf.gradients(x[0], x[1])[0])([loss, I])

    # define model
    return Model([I] + inputs, [grad, loss, feat, bb])


FEAT_SIZE = 100
def setup_single_target(num_models, tokens=[], **kwargs):
    """Vanilla single target model.
    """
    if kwargs.get('clear', False):
        K.clear_session()
        kwargs['clear'] = False
    # loss is -0.25 * log likelihood of observing I_pdb (one-hot encoded)
    # predicted features are positive-definite probabilities (softmax at each pixel)
    def loss_fcn(feat):
        return -0.25 * K.mean(K.sum(I_pdb * K.log(feat + 1e-8), -1), axis=[-1, -2])
    I_pdb = Input(shape=(None, None, FEAT_SIZE), dtype=tf.float32)
    inputs = [I_pdb]
    return setup_model(num_models, loss_fcn, inputs, tokens=tokens, **kwargs)


def setup_dual_target(num_models, tokens=[], **kwargs):
    """Vanilla single target model.
    """
    if kwargs.get('clear', False):
        K.clear_session()
        kwargs['clear'] = False
    def loss(feat):
        loss_1 = -0.25 * K.mean(K.sum(I_pdb_1 * K.log(feat + 1e-8), -1), axis=[-1, -2])
        loss_2 = -0.25 * K.mean(K.sum(I_pdb_2 * K.log(feat + 1e-8), -1), axis=[-1, -2])
        return K.mean([loss_1, loss_2])
    
    I_pdb_1 = Input(shape=(None, None, FEAT_SIZE), dtype=tf.float32)
    I_pdb_2 = Input(shape=(None, None, FEAT_SIZE), dtype=tf.float32)
    inputs = [I_pdb_1, I_pdb_2]
    return setup_model(num_models, loss, inputs, tokens=tokens, **kwargs)


def pssm_to_seq(pssm):
    pssm = np.squeeze(pssm)
    a1 = np.array(alpha_1)
    return ''.join(a1[pssm.argmax(-1)])


def grad_to_df(grad):
    grad = np.squeeze(grad)
    return pd.DataFrame(grad, columns=v12.alpha_1[:-1]).T


def align_pssms(a, b):
    a = pssm_to_seq(a)
    b = pssm_to_seq(b)
    return print_alignment(a, b)


def so_minimizer(opt_rate=2.0, opt_decay=2.0):
    # mutate (apply gradient)
    # allows for negative PSSM
    # averages over all history
    # this is weird.
    def minimizer(pssm, grad, history, opt_iter, **kwargs):
        k = len(history)
        # gradient normalization
        grad_norm = grad / np.linalg.norm(grad, axis=(1, 2), keepdims=True)
        # why not use a tf optimizer like ADAM?
        lr = opt_rate * np.power(1 - k/opt_iter, opt_decay)
        gamma = lr / np.linalg.norm(grad, axis=(1, 2), keepdims=True)
        # print('so_gamma', gamma)
        return pssm - lr * grad_norm
    return minimizer


def mutant_minimizer(num_changes, rule='update all'):
    # more like coordinate descent
    def minimizer(pssm, grad, history, **kwargs):
        try:
            gamma, ix = find_scaling(pssm, -grad, num_changes)
            new_pssm = pssm + -grad*gamma
            print(gamma)
        except ValueError:
            assert False
            return pssm
        if rule == 'update all':
            return new_pssm
        elif rule == 'update one':
            old_seq = pssm[0].argmax(axis=-1)
            new_seq = new_pssm[0].argmax(axis=-1)
            assert (old_seq != new_seq).sum() == num_changes
            print(np.where(old_seq != new_seq)[0])
            return np.eye(20)[new_seq][None]
    return minimizer


def find_scaling(pssm, grad, num_changes, argmax_axis=-1):
    def count_changes(gamma):
        delta = pssm + gamma*grad
        return (pssm.argmax(axis=argmax_axis) != delta.argmax(axis=argmax_axis)).sum()

    def finish(gamma):

        delta = pssm + gamma*grad
        diff = pssm.argmax(axis=argmax_axis) != delta.argmax(axis=argmax_axis)
        diff = np.squeeze(diff)
        ix = np.where(diff)[0]
        return gamma, list(ix)

    guess = 1
    alpha = 2
    stop = 1e-10

    low = guess
    while True:
        changes = count_changes(low)
        if changes == num_changes:
            return finish(low)
        if changes < num_changes:
            break
        low /= alpha

    high = guess * alpha
    while True:
        changes = count_changes(high)
        if changes == num_changes:
            return finish(high)
        if changes > num_changes:
            break
        high *= alpha

    while True:
        test = (high + low)/2
        changes = count_changes(test)
        if changes == num_changes:
            return finish(test)
        if changes < num_changes:
            low = test
        else:
            high = test
            
        if high - low < stop:
            raise ValueError


def design(model, pssm, desired_feat, opt_iter=50, opt_rate=2.0, opt_decay=2.0,
           opt_repeat=1, verbose=False, mask=None):
    # initialize
    if pssm is None:
        pssm = np.random.normal(0, 0.01, size=(1, desired_feat.shape[-3], 20))

    history = []
    for r in range(opt_repeat):
        for k in range(opt_iter):
            # compute gradient
            grad, loss, feat, bb = model.predict([pssm] + desired_feat)
            # gradient normalization
            grad /= np.linalg.norm(grad, axis=(1, 2), keepdims=True)
            # why not use a tf optimizer like ADAM?
            lr = opt_rate * np.power(1 - k/opt_iter, opt_decay)

            # mutate (apply gradient)
            # this allows for negative PSSM
            if mask is not None:
                grad *= mask[None, None]
            pssm -= lr * grad

            if verbose and (k+1) % 5 == 0:
                print(k+1, np.round(loss, 2))

            history += [{'grad': np.squeeze(grad), 
                        'pssm': np.squeeze(pssm), 
                        }]

    return pssm, history


def design2(model, pssm, desired_feats, minimizer, verbose=False, mask=None, 
           opt_iter=50):
    history = []
    # initialize
    if pssm is None:
        pssm = np.random.normal(0, 0.01, size=(1, desired_feat.shape[-3], 20))
    else:
        pssm = pssm.copy()

    for k in range(opt_iter):
        # compute gradient
        grad, loss, feat, bb = model.predict([pssm] + desired_feats)
        if mask is not None:
            grad *= mask[None, None]
        pssm = minimizer(pssm, grad, history, opt_iter=opt_iter)
        
        if verbose and (k+1) % 5 == 0:
            print(k+1, np.round(loss, 2))

        history += [{'grad': np.squeeze(grad), 
                     'pssm': np.squeeze(pssm), 
                    }]
            
    return pssm, history


def dump_histories(histories, filename='misc/sequences.tif', grad_cutoff=0):
    import sys
    sys.path.append('packages/OpticalPooledScreens/')
    from ops.io import save_stack, GLASBEY
    arr = []
    for history in histories:
        for h in history:
            A = h['pssm'].argmax(-1)
            B = -h['grad'].copy()
            B[B < grad_cutoff] = 0
            B = B.argmax(-1)
            arr += [[A, B]]

    s = np.array(arr).transpose([1, 0, 2])
    save_stack(filename, s, luts=(GLASBEY, GLASBEY))