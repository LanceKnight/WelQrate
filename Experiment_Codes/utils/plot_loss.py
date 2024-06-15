import matplotlib.pyplot as plt
import os

def plot_epoch():
    # Define the filename and load the data
    filename = "logs/loss_per_epoch.log"
    with open(filename, "r") as f:
        lines = f.readlines()

    epochs = []
    losses = []
    logAUCs = []
    lrs =[]
    # ts = []

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 5))
    # Parse the data from each line
    for line in lines:

        # epoch, loss,  lr, t = line.strip().split("\t")
        if len(line.strip().split("\t"))==3:
            epoch, loss, lr = line.strip().split("\t")
        elif len(line.strip().split("\t"))==4:
            epoch, loss, logAUC, lr = line.strip().split("\t")
            logAUCs.append(float(logAUC.split("=")[1]))
        epochs.append(int(epoch))
        losses.append(float(loss.split("=")[1]))
        lrs.append(float(lr.split("=")[1]))
        # ts.append(float(t.split("=")[1]))

    # Plot the loss vs epoch
    ax1.scatter(epochs, losses)
    ax1.set_xlabel("Epoch")
    ax1.set_ylabel("Loss")

    # Plot the logAUC vs epoch
    if len(line.strip().split("\t")) == 4:
        print(f'epoch={epoch}, logAUC={logAUCs}')
        ax2.scatter(epochs, logAUCs)
    ax2.set_xlabel("Epoch")
    ax2.set_ylabel("logAUC")

    ax3.scatter(epochs, lrs)
    ax3.set_xlabel("Epoch")
    ax3.set_ylabel("learning rate")

    # ax4.scatter(epochs, ts)
    # ax4.set_xlabel("Epoch")
    # ax4.set_ylabel("temperature")

    # Save the plot as a PNG file
    os.makedirs('images', exist_ok=True)
    plt.savefig("images/loss_and_auc_vs_epoch.png")

if __name__ == "__main__":
    plot_epoch()