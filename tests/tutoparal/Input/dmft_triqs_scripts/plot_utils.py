import matplotlib.pyplot as plt
import tkinter as tk

HA_EV = 27.211386245981

class Figure:

    def __init__(self):

        # --- Detect screen dimensions and DPI ---
        root = tk.Tk()
        root.withdraw()
        screen_w = root.winfo_screenwidth()
        screen_h = root.winfo_screenheight()
        dpi      = root.winfo_fpixels('1i')
        root.destroy()

        self.screen_w, self.screen_h = screen_w, screen_h

        # --- Compute scaling factors ---
        font_ratio   = screen_w * dpi / (3440 * 96)
        aspect_ratio = 1093 / 2400
        scale_ratio  = min(1.0, screen_h / (screen_w * aspect_ratio))

        self.font_ratio = font_ratio

        target_w = int(screen_w * scale_ratio)
        target_h = int(target_w * aspect_ratio)

        self.target_w, self.target_h = target_w, target_h

        fig_width  = target_w / dpi
        fig_height = target_h / dpi

        # --- Plot settings ---
        plt.rcParams.update({
            "axes.edgecolor": "black",
            "axes.linewidth": 5 * font_ratio
        })

        fig, ax = plt.subplots(figsize=(fig_width, fig_height), dpi=dpi)

        self.ax = ax

        ax.tick_params(
            axis="both",
            which="major",
            labelsize=50 * font_ratio,
            width=5 * font_ratio,
            length=20 * font_ratio,
            pad=15 * font_ratio
        )

        ax.yaxis.set_label_coords(-0.12, 0.45)
        ax.xaxis.set_label_coords(0.55, -0.2)

        fig.subplots_adjust(bottom=0.22, top=0.95, left=0.15, right=0.9)

    def add_data(self, x, y, params):

        ax = self.ax
        font_ratio = self.font_ratio

        if "yerr" in params:
            ax.errorbar(x, y, params["yerr"], lw=2*font_ratio, color="black", ecolor="#c1272d", elinewidth=5*font_ratio)
        else:
            ax.plot(x, y, lw=5*font_ratio, color="#c1272d")

        if "xlabel" in params:
            ax.set_xlabel(params["xlabel"], fontsize=50*font_ratio)

        if "ylabel" in params:
            ax.set_ylabel(params["ylabel"], fontsize=50*font_ratio)

        if "xlim" in params:
            ax.set_xlim(params["xlim"][0],params["xlim"][1])

        if "xticks" in params:
            ax.set_xticks(params["xticks"])

        if "xticklabels" in params:
            ax.set_xticklabels(params["xticklabels"])

    def plot(self):

        target_w, target_h = self.target_w, self.target_h
        screen_w, screen_h = self.screen_w, self.screen_h

        # --- Actually resize or maximize the window ---
        try:
            manager = plt.get_current_fig_manager()
            backend = plt.get_backend().lower()

            if "qt" in backend:
                manager.window.showNormal()
                manager.window.resize(target_w, target_h)
                manager.window.move((screen_w - target_w) // 2, (screen_h - target_h) // 2)

            elif "tk" in backend:
                manager.window.wm_geometry(f"{target_w}x{target_h}+{(screen_w - target_w)//2}+{(screen_h - target_h)//2}")

            elif "wx" in backend:
                manager.window.SetSize((target_w, target_h))

        except Exception as e:
            print(f"[Warning] Could not resize window automatically: {e}")

        plt.show()

