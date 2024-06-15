import numpy as np
from sklearn.metrics import confusion_matrix, auc, roc_curve, f1_score, \
    roc_auc_score


def sigmoid(x):
    z = 1 / (1 + np.exp(-x))
    return z

def MAE(true_y, predicted_score):
    return np.mean(np.abs(true_y - predicted_score))

def MSE(true_y, predicted_score):
    return np.mean((true_y - predicted_score) ** 2)

def RMSE(true_y, predicted_score):
    return np.sqrt(np.mean((true_y - predicted_score) ** 2))

def R2(true_y, predicted_score):
    SS_res = np.sum((true_y - predicted_score) ** 2)
    SS_tot = np.sum((true_y - np.mean(true_y)) ** 2)
    return 1 - SS_res / SS_tot



def calculate_logAUC(true_y, predicted_score, FPR_range=(0.001, 0.1)):
    """
    Author: Yunchao "Lance" Liu (lanceknight26@gmail.com)
    Calculate logAUC in a certain FPR range (default range: [0.001, 0.1]).
    This was used by previous methods [1] and the reason is that only a
    small percentage of samples can be selected for experimental tests in
    consideration of cost. This means only molecules with very high
    predicted activity values can be worth testing, i.e., the decision
    threshold is high. And the high decision threshold corresponds to the
    left side of the ROC curve, i.e., those FPRs with small values. Also,
    because the threshold cannot be predetermined, the area under the curve
    is used to consolidate all possible thresholds within a certain FPR
    range. Finally, the logarithm is used to bias smaller FPRs. The higher
    the logAUC[0.001, 0.1], the better the performance.

    A perfect classifer gets a logAUC[0.001, 0.1] ) of 1, while a random
    classifer gets a logAUC[0.001, 0.1] ) of around 0.0215 (See [2])

    References:
    [1] Mysinger, M.M. and B.K. Shoichet, Rapid Context-Dependent Ligand
    Desolvation in Molecular Docking. Journal of Chemical Information and
    Modeling, 2010. 50(9): p. 1561-1573.
    [2] Mendenhall, J. and J. Meiler, Improving quantitative
    structure–activity relationship models using Artificial Neural Networks
    trained with dropout. Journal of computer-aided molecular design,
    2016. 30(2): p. 177-189.
    :param true_y: numpy array of the ground truth. Values are either 0 (
    inactive) or 1(active).
    :param predicted_score: numpy array of the predicted score (The
    score does not have to be between 0 and 1)
    :param FPR_range: the range for calculating the logAUC formated in
    (x, y) with x being the lower bound and y being the upper bound
    :return: a numpy array of logAUC of size [1,1]
    """

    np.seterr(divide='ignore')
    # FPR range validity check
    if FPR_range == None:
        raise Exception('FPR range cannot be None')
    lower_bound = FPR_range[0]
    upper_bound = FPR_range[1]
    if (lower_bound >= upper_bound):
        raise Exception('FPR upper_bound must be greater than lower_bound')

    fpr, tpr, thresholds = roc_curve(true_y, predicted_score, pos_label=1)

    tpr = np.append(tpr, np.interp([lower_bound, upper_bound], fpr, tpr))
    fpr = np.append(fpr, [lower_bound, upper_bound])

    # Sort both x-, y-coordinates array
    tpr = np.sort(tpr)
    fpr = np.sort(fpr)

    # Get the data points' coordinates. log_fpr is the x coordinate, tpr is the y coordinate.
    log_fpr = np.log10(fpr)
    x = log_fpr
    y = tpr
    lower_bound = np.log10(lower_bound)
    upper_bound = np.log10(upper_bound)

    # Get the index of the lower and upper bounds
    lower_bound_idx = np.where(x == lower_bound)[-1][-1]
    upper_bound_idx = np.where(x == upper_bound)[-1][-1]

    # Create a new array trimmed at the lower and upper bound
    trim_x = x[lower_bound_idx:upper_bound_idx + 1]
    trim_y = y[lower_bound_idx:upper_bound_idx + 1]

    area = auc(trim_x, trim_y) / (upper_bound - lower_bound)
    return area



def cal_EF(true_y, predicted_score, k):
    # Ns: the number of compounds in the selection set (the predicted actives), k
    # N: the number of compounds in the entire dataset
    # ns: the number of true active compounds in the selection set
    # n: the number of true active compounds in the entire dataset

    order_idx_pred = np.argsort(-predicted_score)

    order_predicted_score = predicted_score[order_idx_pred]
    order_true_y = true_y[order_idx_pred]

    N = len(true_y)
    Ns = k

    ns = (order_true_y[:k] == 1).sum()
    n = true_y.sum()

    return N * ns / (Ns * n)



def cal_DCG(true_y, predicted_score, k):
    # assuming k<= len(true_y) == len(predicted_score)

    order_idx_pred = np.argsort(-predicted_score)
    order_idx_label = np.argsort(-true_y)

    dcg = np.sum(true_y[order_idx_pred][:k] / np.log2(np.arange(1 + 1, (k + 1) + 1)))

    return dcg


def cal_BEDROC_score(true_y, predicted_score, decreasing=True, alpha=20.0):
    """BEDROC metric implemented according to Truchon and Bayley.

    The Boltzmann Enhanced Descrimination of the Receiver Operator
    Characteristic (BEDROC) score is a modification of the Receiver Operator
    Characteristic (ROC) score that allows for a factor of *early recognition*.

    References:
        The original paper by Truchon et al. is located at `10.1021/ci600426e
        <http://dx.doi.org/10.1021/ci600426e>`_.

    Args:
        true_y (array_like):
            Binary class labels. 1 for positive class, 0 otherwise.
        predicted_score (array_like):
            Prediction values.
        decreasing (bool):
            True if high values of ``predicted_score`` correlates to positive class.
        alpha (float):
            Early recognition parameter.

    Returns:
        float:
            Value in interval [0, 1] indicating degree to which the predictive
            technique employed detects (early) the positive class.
     """

    assert len(true_y) == len(predicted_score), \
        'The number of scores must be equal to the number of labels'

    big_n = len(true_y)
    n = sum(true_y == 1)

    if decreasing:
        order = np.argsort(-predicted_score)
    else:
        order = np.argsort(predicted_score)

    m_rank = (true_y[order] == 1).nonzero()[0]

    s = np.sum(np.exp(-alpha * m_rank / big_n))

    r_a = n / big_n

    rand_sum = r_a * (1 - np.exp(-alpha)) / (np.exp(alpha / big_n) - 1)

    fac = r_a * np.sinh(alpha / 2) / (np.cosh(alpha / 2) -
                                      np.cosh(alpha / 2 - alpha * r_a))

    cte = 1 / (1 - np.exp(alpha * (1 - r_a)))

    return s * fac / rand_sum + cte