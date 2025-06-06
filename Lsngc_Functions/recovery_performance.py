def recovery_performance(Adjs,label):
    Adjs=Adjs.copy()
    label=label.copy()
    N=len(Adjs)
    auc_all = np.zeros((N), dtype=np.float32)   
    for i in range(N):
        auc_all[i] = metrics.roc_auc_score(label[i].flatten(), Adjs[i].flatten())
    return auc_all
