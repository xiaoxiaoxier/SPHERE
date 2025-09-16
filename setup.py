from setuptools import setup, find_packages

setup(
    name="sphere",  # 包名（pip install sphere）
    version="0.1.0",
    author="Your Name",
    author_email="your.email@example.com",
    description="SPHERE: Spatial Proximity-Weighted High-resolution Expression Reconstruction for Visium HD",
    long_description=open("README.md", encoding="utf-8").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/<your-username>/SPHERE",  # 改成你的 GitHub 仓库地址
    license="MIT",
    packages=find_packages(),  # 自动找到 sphere/ 这个包
    include_package_data=True,
    install_requires=open("requirements.txt", encoding="utf-8").read().splitlines(),
    python_requires=">=3.8",  # 根据你的环境适当调整
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)