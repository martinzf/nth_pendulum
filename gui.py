import numpy as np
import tkinter as tk
import customtkinter as ctk

NMAX = 50
WIDTH, HEIGHT = 500, 450 # pixels

class App(ctk.CTk):
    def __init__(self):
        super().__init__()
        # Styling
        self.title('Simulation parameters')    
        self.geometry(f'{WIDTH}x{HEIGHT}') 
        self.frame = ctk.CTkFrame(self) 
        self.frame.pack(padx=20, pady=20, fill='both', expand=True)
        self.frame.columnconfigure(0, weight=1)
        self.frame.columnconfigure(1, weight=2)
        self.frame.columnconfigure(2, weight=2)

        # n
        numberLabel = ctk.CTkLabel(self.frame, text='n = 1')
        numberLabel.grid(row=0, column=0, padx=10, pady=10, sticky='ew')
        number = tk.IntVar(value=1)
        def update_label(_):
            numberLabel.configure(text=f'n = {number.get()}')
        numberSlider = ctk.CTkSlider(self.frame, from_=1, to=NMAX, variable=number, command=update_label)
        numberSlider.grid(row=0, column=1, columnspan=2, padx=10, pady=10, sticky='ew')

        # theta0 
        self.theta0Label = ctk.CTkLabel(self.frame, text='θₒ (rad)')
        self.theta0Label.grid(row=1, column=0, padx=10, pady=10, sticky='ew')
        self.theta0Entry = ctk.CTkEntry(self.frame, 
                                   placeholder_text='Expected real NumPy (n,) ndarray', 
                                   placeholder_text_color='red')
        self.theta0Entry.insert(0, 'np.ones(n)')
        self.theta0Entry.grid(row=1, column=1, columnspan=3, padx=10, pady=10, sticky='ew')

        # d_theta0 
        self.d_theta0Label = ctk.CTkLabel(self.frame, text='dθₒ/dt (rad/s)')
        self.d_theta0Label.grid(row=2, column=0, padx=10, pady=10, sticky='ew')
        self.d_theta0Entry = ctk.CTkEntry(self.frame, 
                                     placeholder_text='Expected real NumPy (n,) ndarray', 
                                     placeholder_text_color='red')
        self.d_theta0Entry.insert(0, 'np.zeros(n)')
        self.d_theta0Entry.grid(row=2, column=1, columnspan=3, padx=10, pady=10, sticky='ew')

        # l 
        self.lLabel = ctk.CTkLabel(self.frame, text='ℓ (m)')
        self.lLabel.grid(row=3, column=0, padx=10, pady=10, sticky='ew')
        self.lEntry = ctk.CTkEntry(self.frame, 
                                   placeholder_text='Expected strictly positive NumPy (n,) ndarray', 
                                   placeholder_text_color='red')
        self.lEntry.insert(0, 'np.ones(n)')
        self.lEntry.grid(row=3, column=1, columnspan=3, padx=10, pady=10, sticky='ew')

        # m 
        self.mLabel = ctk.CTkLabel(self.frame, text='m (kg)')
        self.mLabel.grid(row=4, column=0, padx=10, pady=10, sticky='ew')
        self.mEntry = ctk.CTkEntry(self.frame, 
                                 placeholder_text='Expected strictly positive NumPy (n,) ndarray', 
                                 placeholder_text_color='red')
        self.mEntry.insert(0, 'np.ones(n)')
        self.mEntry.grid(row=4, column=1, columnspan=3, padx=10, pady=10, sticky='ew')

        # g 
        self.gLabel = ctk.CTkLabel(self.frame, text='g (m/s²)')
        self.gLabel.grid(row=5, column=0, padx=10, pady=10, sticky='ew')
        self.gEntry = ctk.CTkEntry(self.frame, 
                              placeholder_text='Expected int or float', 
                              placeholder_text_color='red')
        self.gEntry.insert(0, '9.8')
        self.gEntry.grid(row=5, column=1, columnspan=3, padx=10, pady=10, sticky='ew')

        # Duration
        self.durLabel = ctk.CTkLabel(self.frame, text='t (s)')
        self.durLabel.grid(row=6, column=0, padx=10, pady=10, sticky='ew')
        self.durEntry = ctk.CTkEntry(self.frame, 
                              placeholder_text='Expected strictly positive int or float', 
                              placeholder_text_color='red')
        self.durEntry.insert(0, '20')
        self.durEntry.grid(row=6, column=1, columnspan=3, padx=10, pady=10, sticky='ew')

        # Simulate Button
        def get_data():
            global n, theta0, d_theta0, l, m, g, dur
            n = number.get()
            theta0 = check_vector(n, self.theta0Entry)
            d_theta0 = check_vector(n, self.d_theta0Entry)
            l = check_vector(n, self.lEntry)
            m = check_vector(n, self.mEntry)
            g = check_scalar(self.gEntry)
            dur = check_scalar(self.durEntry)
            if None in [*theta0, *d_theta0, *l, *m, g, dur]:
                return
            if np.any(l <= 0):
                self.lEntry.delete(0, 'end')
                self.focus()
                return
            if np.any(m <= 0):
                self.mEntry.delete(0, 'end')
                self.focus()
                return
            if dur <= 0:
                self.durEntry.delete(0, 'end')
                self.focus()
                return
            self.quit()

        def check_vector(n, xEntry):
            try:
                x = eval(xEntry.get())
                if not(isinstance(x, np.ndarray) and np.shape(x)!=(n)) or x.dtype == complex:
                    xEntry.delete(0, 'end')
                    self.focus()
                    return 
                return x
            except:
                xEntry.delete(0, 'end')
                self.focus()
                return 

        def check_scalar(xEntry):
            try:
                x = eval(xEntry.get())
                if not(isinstance(x, (int, float))):
                    xEntry.delete(0, 'end')
                    self.focus()
                    return 
                return x
            except:
                xEntry.delete(0, 'end')
                self.focus()
                return 

        self.SimulateButton = ctk.CTkButton(self.frame, text='Simulate', command=get_data)
        self.SimulateButton.grid(row=7, column=1, padx=20, pady=20, sticky='ew')

def setup():
    ctk.set_appearance_mode("System")  
    ctk.set_default_color_theme("blue") 
    app = App()
    app.mainloop()
    app.quit()
    return n, theta0, d_theta0, l, m, g, dur