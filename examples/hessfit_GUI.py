import tkinter as tk
from tkinter import filedialog, messagebox
import subprocess
import os

def print_init():
    info = """
 ======================================================
   Program:      hessfit
   Creator:      Emanuele Falbo, Napoli
   Language:     Python 3.v later
   Description:  The program returns force constants for
                 bonded values and non-bonded parameters
   Mail:         falbo.emanuele@gmail.com 
 =======================================================
    """
    print(info)
    messagebox.showinfo("Program Info", info)

def run_hessfit(optfile, g09path):
    try:
        if not g09path or not os.path.exists(os.path.join(g09path, "g09")):
            messagebox.showerror("Error", "Invalid Gaussian path")
            return

        SM = "hessfit_harmonic.py"
        BS = "build_4_hessfit.py"

        # Run build script
        subprocess.run([BS, optfile], check=True)

        # Run Gaussian on required files
        for f in ["GauHarm.gjf", "GauNonBon.gjf"]:
            subprocess.run([os.path.join(g09path, "g09"), f], check=True)

        # Format check the .chk files
        for f in ["GauHarm.chk", "GauNonBon.chk"]:
            subprocess.run([os.path.join(g09path, "formchk"), "-3", f, f"{os.path.splitext(f)[0]}.fchk"], check=True)

        # Execute Harmonic script
        subprocess.run([SM, optfile], check=True)

        # Final Gaussian execution
        subprocess.run([os.path.join(g09path, "g09"), "hessfit4gau.gjf"], check=True)

        messagebox.showinfo("Success", "hessfit execution completed!")
    except subprocess.CalledProcessError as e:
        messagebox.showerror("Error", f"Subprocess failed: {e}")
    except Exception as e:
        messagebox.showerror("Error", f"An error occurred: {e}")

def browse_file(entry):
    file_path = filedialog.askopenfilename()
    if file_path:
        entry.delete(0, tk.END)
        entry.insert(0, file_path)

def browse_directory(entry):
    dir_path = filedialog.askdirectory()
    if dir_path:
        entry.delete(0, tk.END)
        entry.insert(0, dir_path)

def main():
    print_init()

    # Create main window
    root = tk.Tk()
    root.title("hessfit GUI")

    # Optfile input
    tk.Label(root, text="Optfile (.json):").grid(row=0, column=0, padx=10, pady=5)
    optfile_entry = tk.Entry(root, width=50)
    optfile_entry.grid(row=0, column=1, padx=10, pady=5)
    tk.Button(root, text="Browse", command=lambda: browse_file(optfile_entry)).grid(row=0, column=2, padx=10, pady=5)

    # Gaussian path input
    tk.Label(root, text="Gaussian Path:").grid(row=1, column=0, padx=10, pady=5)
    g09path_entry = tk.Entry(root, width=50)
    g09path_entry.grid(row=1, column=1, padx=10, pady=5)
    tk.Button(root, text="Browse", command=lambda: browse_directory(g09path_entry)).grid(row=1, column=2, padx=10, pady=5)

    # Run button
    tk.Button(root, text="Run hessfit", command=lambda: run_hessfit(optfile_entry.get(), g09path_entry.get())).grid(
        row=2, column=1, pady=20
    )

    # Start GUI loop
    root.mainloop()

if __name__ == "__main__":
    main()
