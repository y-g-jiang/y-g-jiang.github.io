import tkinter as tk
from tkinter import filedialog, messagebox, ttk
from PIL import Image, ImageTk
import math
import cv2
import numpy as np
import webbrowser  # 导入 webbrowser 模块

class ImageViewer(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("移轴畸变及画面边角畸变修正器v1.0.1 by姜尧耕")
        # self.wm_state('zoomed')  # 可以移除或注释掉

        # 获取屏幕宽度和高度
        screen_width = self.winfo_screenwidth()
        screen_height = self.winfo_screenheight()
        self.geometry(f"{screen_width}x{screen_height}+0+0") # 设置窗口大小和位置

        # 移除窗口图标 (跨平台方法)
        self.withdraw()  # 隐藏主窗口
        self.iconphoto(False, tk.PhotoImage(width=1, height=1))
        self.deiconify() # 重新显示窗口

        # 使用 ttk 样式
        self.style = ttk.Style(self)
        self.style.theme_use('clam')

        # 设置默认字体
        self.default_font = ('Helvetica', 10)
        self.style.configure('.', font=self.default_font, background='#f0f0f0')
        self.style.configure('TButton', padding=5)
        self.style.configure('TLabel', background='#f0f0f0')
        self.style.configure('TFrame', background='#f0f0f0')

        self.original_image = None
        self.image = None
        self.tk_image = None

        # 主内容 Frame
        self.content_frame = ttk.Frame(self, padding=(10, 10, 10, 10))
        self.content_frame.pack(fill=tk.BOTH, expand=True)

        # 左侧显示图片的 Frame
        self.image_frame = ttk.Frame(self.content_frame)
        self.image_frame.grid(row=0, column=0, sticky='nsew', padx=(0, 10))
        self.image_frame.bind("<Configure>", self.on_frame_resize)

        self.image_label = ttk.Label(self.image_frame)
        self.image_label.pack(anchor="center", fill=tk.BOTH, expand=True)
        self.image_label.bind("<Button-1>", self.on_image_click)

        # 右侧输入框和按钮的 Frame
        self.button_frame = ttk.Frame(self.content_frame)
        self.button_frame.grid(row=0, column=1, sticky='ns')

        # 等效全幅焦距
        focal_length_frame = ttk.Frame(self.button_frame)
        focal_length_frame.pack(pady=(10, 5), fill=tk.X)
        self.focal_length_label = ttk.Label(focal_length_frame, text="等效全幅焦距 (mm)")
        self.focal_length_label.pack(side=tk.LEFT)
        self.focal_length_entry = ttk.Entry(focal_length_frame)
        self.focal_length_entry.pack(side=tk.RIGHT, expand=True, fill=tk.X)

        # 横竖比例
        aspect_ratio_frame = ttk.Frame(self.button_frame)
        aspect_ratio_frame.pack(pady=(5, 10), fill=tk.X)
        self.aspect_ratio_label = ttk.Label(aspect_ratio_frame, text="横竖比例")
        self.aspect_ratio_label.pack(side=tk.LEFT)
        self.bi1_entry = ttk.Entry(aspect_ratio_frame, width=5, justify='center')
        self.bi1_entry.pack(side=tk.LEFT)
        self.aspect_ratio_separator = ttk.Label(aspect_ratio_frame, text=":")
        self.aspect_ratio_separator.pack(side=tk.LEFT)
        self.bi2_entry = ttk.Entry(aspect_ratio_frame, width=5, justify='center')
        self.bi2_entry.pack(side=tk.LEFT)

        # 导入图片按钮
        self.load_button = ttk.Button(self.button_frame, text="导入图片", command=self.load_image)
        self.load_button.pack(pady=(10, 10), fill=tk.X)

        self.coordinate_label = ttk.Label(self.button_frame, text="点击图片显示坐标比例")
        self.coordinate_label.pack(pady=(5, 10), fill=tk.X)

        self.apply_button = ttk.Button(self.button_frame, text="应用透视变换", command=self.apply_perspective_transform, state=tk.DISABLED)
        self.apply_button.pack(pady=(10, 10), fill=tk.X)

        # 重置图片按钮
        self.reset_button = ttk.Button(self.button_frame, text="完全重置图片", command=self.reset_image, state=tk.DISABLED)
        self.reset_button.pack(pady=(10, 10), fill=tk.X)

        # 导出图片按钮
        self.export_button = ttk.Button(self.button_frame, text="导出清晰图片", command=self.export_image, state=tk.DISABLED)
        self.export_button.pack(pady=(10, 10), fill=tk.X)

        # 版权信息
        self.copyright_label = ttk.Label(self.button_frame, text="by 姜尧耕，网上搜我名字就行")
        self.copyright_label.pack(pady=(10, 5), fill=tk.X)

        # 备用信息
        self.fallback_label = ttk.Label(self.button_frame, text="如果链接失效，网上搜索一下软件名`移轴畸变及画面边角畸变修正器`即可。")
        self.fallback_label.pack(pady=(0, 5), fill=tk.X)

        # 打开百度按钮
        self.baidu_button = ttk.Button(self.button_frame, text="使用教程及版本更新", command=self.open_baidu)
        self.baidu_button.pack(pady=(5, 10), fill=tk.X)

        # 配置行和列的权重，使图片 Frame 可以扩展
        self.content_frame.columnconfigure(0, weight=1)
        self.content_frame.rowconfigure(0, weight=1)

        self.biliheng = None
        self.bilishu = None

    def open_baidu(self):
        webbrowser.open("https://zhuanlan.zhihu.com/p/18525694310")

    def load_image(self):
        file_path = filedialog.askopenfilename(filetypes=[("Image files", "*.png;*.jpg;*.jpeg;*.gif;*.bmp")])
        if file_path:
            try:
                self.original_image = Image.open(file_path)
                self.image = self.original_image.copy() # 保存原始图像
                self.display_image()
                self.apply_button.config(state=tk.DISABLED)
                self.reset_button.config(state=tk.DISABLED)
                self.export_button.config(state=tk.NORMAL)
                self.coordinate_label.config(text="点击图片显示坐标比例")
                self.biliheng = None
                self.bilishu = None
            except Exception as e:
                print(f"Error loading image: {e}")
                messagebox.showerror("错误", f"无法加载图片: {e}")

    def display_image(self):
        if self.image:
            frame_width = self.image_frame.winfo_width()
            frame_height = self.image_frame.winfo_height()
            img_width, img_height = self.image.size

            # 计算在不改变图片比例的情况下，适应 frame 的缩放比例
            width_ratio = frame_width / img_width
            height_ratio = frame_height / img_height
            scale = min(width_ratio, height_ratio)

            # 如果图片比 frame 大，则缩小以适应
            if scale < 1:
                new_width = int(img_width * scale)
                new_height = int(img_height * scale)
                resized_image = self.image.resize((new_width, new_height), Image.LANCZOS)
            else:
                # 否则，保持原始分辨率
                resized_image = self.image

            self.tk_image = ImageTk.PhotoImage(resized_image)
            self.image_label.config(image=self.tk_image)
            self.image_label.image = self.tk_image

    def on_frame_resize(self, event):
        if self.image:
            self.display_image()

    def on_image_click(self, event):
        if self.image:
            displayed_width = self.image_label.winfo_width()
            displayed_height = self.image_label.winfo_height()

            if displayed_width == 0 or displayed_height == 0:
                return  # 避免除以零

            center_x = displayed_width / 2
            center_y = displayed_height / 2

            offset_x = event.x - center_x
            offset_y = event.y - center_y

            self.biliheng = offset_x / center_x
            self.bilishu = -offset_y / center_y

            self.coordinate_label.config(text=f"X比例: {self.biliheng:.4f}, Y比例: {self.bilishu:.4f}")
            if self.original_image:
                self.apply_button.config(state=tk.NORMAL)
                self.reset_button.config(state=tk.NORMAL)

    def apply_perspective_transform(self):
        if self.original_image and self.biliheng is not None and self.bilishu is not None:
            try:
                bi1 = float(self.bi1_entry.get())
                bi2 = float(self.bi2_entry.get())
                equivalent_focal_length = float(self.focal_length_entry.get())

                if bi1 <= 0 or bi2 <= 0 or equivalent_focal_length <= 0:
                    messagebox.showerror("错误", "横竖比例和等效全幅焦距必须是正数。")
                    return

                halfd = 1 * math.tan(math.atan(43.4 / (2 * equivalent_focal_length)))

                ax = (bi1 / math.sqrt(bi1**2 + bi2**2)) * halfd
                ay = (bi2 / math.sqrt(bi1**2 + bi2**2)) * halfd

                def Mperspectivetransform(wA, hA, sub_a_x, sub_a_y):
                    numerator1 = math.sqrt(1 + hA**2 + wA**2) * (-wA + sub_a_x)
                    denominator1 = math.sqrt(1 + wA**2) * (1 + wA * sub_a_x + hA * sub_a_y)
                    term1 = numerator1 / denominator1

                    numerator2 = -hA - hA * wA * sub_a_x + (1 + wA**2) * sub_a_y
                    denominator2 = math.sqrt(1 + wA**2) * (1 + wA * sub_a_x + hA * sub_a_y)
                    term2 = numerator2 / denominator2

                    return [term1, term2]

                M1_values = {'sub_a_x': ax, 'sub_a_y': ay}
                M2_values = {'sub_a_x': -ax, 'sub_a_y': ay}
                M3_values = {'sub_a_x': ax, 'sub_a_y': -ay}
                M4_values = {'sub_a_x': -ax, 'sub_a_y': -ay}

                wA_sub = self.biliheng * ax
                hA_sub = self.bilishu * ay

                M1 = Mperspectivetransform(wA_sub, hA_sub, M1_values['sub_a_x'], M1_values['sub_a_y'])
                M2 = Mperspectivetransform(wA_sub, hA_sub, M2_values['sub_a_x'], M2_values['sub_a_y'])
                M3 = Mperspectivetransform(wA_sub, hA_sub, M3_values['sub_a_x'], M3_values['sub_a_y'])
                M4 = Mperspectivetransform(wA_sub, hA_sub, M4_values['sub_a_x'], M4_values['sub_a_y'])

                img_width, img_height = self.original_image.size
                source_points = np.float32([[0, 0], [img_width, 0], [img_width, img_height], [0, img_height]])

                # 注意这里的顺序：M2(左上), M1(右上), M3(右下), M4(左下)
                dest_points = np.float32([
                    [(M2[0] + 1) * img_width / 2, (1 - M2[1]) * img_height / 2],  # 左上
                    [(M1[0] + 1) * img_width / 2, (1 - M1[1]) * img_height / 2],  # 右上
                    [(M3[0] + 1) * img_width / 2, (1 - M3[1]) * img_height / 2],  # 右下
                    [(M4[0] + 1) * img_width / 2, (1 - M4[1]) * img_height / 2]   # 左下
                ])

                M = cv2.getPerspectiveTransform(source_points, dest_points)
                warped_image_np = cv2.warpPerspective(np.array(self.original_image), M, (img_width, img_height))
                transformed_image = Image.fromarray(warped_image_np)

                # 获取透视变换后图像的尺寸
                transformed_width, transformed_height = transformed_image.size

                # 应用比例压缩
                bi1 = float(self.bi1_entry.get())
                bi2 = float(self.bi2_entry.get())
                compression_ratio = bi1 / bi2
                new_height = int(transformed_height * compression_ratio)
                self.image = transformed_image.resize((transformed_width, new_height), Image.LANCZOS)

                self.display_image()
                self.apply_button.config(state=tk.DISABLED)
                # 在应用透视变换后启用重置按钮
                self.reset_button.config(state=tk.NORMAL)
                self.export_button.config(state=tk.NORMAL)

            except ValueError:
                messagebox.showerror("错误", "请输入有效的数字。")
            except Exception as e:
                messagebox.showerror("错误", f"应用透视变换时出错: {e}")

    def reset_image(self):
        if self.original_image:
            self.image = self.original_image.copy()
            self.display_image()
            self.apply_button.config(state=tk.DISABLED)
            self.reset_button.config(state=tk.DISABLED)
            self.export_button.config(state=tk.NORMAL)
            self.coordinate_label.config(text="点击图片显示坐标比例")
            self.biliheng = None
            self.bilishu = None

    def export_image(self):
        if self.image:
            file_path = filedialog.asksaveasfilename(defaultextension=".png",
                                                     filetypes=[("PNG files", "*.png"),
                                                                ("JPEG files", "*.jpg"),
                                                                ("BMP files", "*.bmp"),
                                                                ("All files", "*.*")])
            if file_path:
                try:
                    self.image.save(file_path)
                    messagebox.showinfo("导出成功", f"图片已成功导出到: {file_path}")
                except Exception as e:
                    messagebox.showerror("导出失败", f"导出图片时出错: {e}")
        else:
            messagebox.showerror("错误", "没有可导出的图片。请先加载图片。")

if __name__ == "__main__":
    app = ImageViewer()
    app.mainloop()
