# Adapted from the scheduler used in Graphormer:
# https://github.com/microsoft/Graphormer/tree/ogb-lsc/OGB-LSC/graphormer/src
# Copyright (c) Microsoft Corporation.
# Licensed under the MIT License.

import math
from torch.optim.lr_scheduler import _LRScheduler


class PolynomialDecayLR(_LRScheduler):

    def __init__(self, optimizer, warmup_iterations, tot_iterations, lr, end_lr, power, last_epoch=-1, verbose=False):
        self.warmup_updates = warmup_iterations
        self.tot_updates = tot_iterations
        self.lr = lr
        self.end_lr = end_lr
        self.power = power
        super(PolynomialDecayLR, self).__init__(optimizer, last_epoch, verbose)

    def get_lr(self):
        if self._step_count <= self.warmup_updates:
            self.warmup_factor = self._step_count / float(self.warmup_updates)
            lr = self.warmup_factor * self.lr
        elif self._step_count >= self.tot_updates:
            lr = self.end_lr
        else:
            warmup = self.warmup_updates
            lr_range = self.lr - self.end_lr
            pct_remaining = 1 - (self._step_count - warmup) / (
                self.tot_updates - warmup
            )
            lr = lr_range * pct_remaining ** (self.power) + self.end_lr

        return [lr for group in self.optimizer.param_groups]

    def _get_closed_form_lr(self):
        assert False

def get_scheduler(optimizer, config, dataset):

    num_train_batches = math.ceil(len(dataset)/int(config['TRAIN']['batch_size']))

    tot_iterations = (num_train_batches) * int(config['TRAIN']['num_epochs'])
    print(f'tot_iterations={tot_iterations}')

    scheduler = {
        'scheduler': PolynomialDecayLR(
            optimizer,
            warmup_iterations=int(config['TRAIN']['warmup_iterations']),
            tot_iterations=tot_iterations,
            lr=float(config['TRAIN']['peak_lr']),
            end_lr=float(config['TRAIN']['end_lr']),
            power=1.0,
        ),
        'dataset_name': 'learning_rate',
        'interval': 'step',
        'frequency': 1,
    }
    return scheduler['scheduler']

def get_lr(optimizer):
    for param_group in optimizer.param_groups:
        return param_group['lr']