from selenium import webdriver
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.chrome.options import Options
from webdriver_manager.chrome import ChromeDriverManager
from selenium.webdriver.common.by import By
import matplotlib.pyplot as plt
from PIL import Image
import requests
from io import BytesIO
import time

# Set Chrome options
chrome_options = Options()
chrome_options.add_argument("--headless")  # Run in headless mode

# Automatically use the right ChromeDriver
service = Service(ChromeDriverManager().install())
driver = webdriver.Chrome(service=service, options=chrome_options)

# Open NASA Earth Observatory
driver.get("https://earthobservatory.nasa.gov/images")

# Wait and grab image
time.sleep(2)
first_image = driver.find_element(By.CSS_SELECTOR, "div.item a img")
img_url = first_image.get_attribute("src")

driver.quit()

# Download and display image
response = requests.get(img_url)
img = Image.open(BytesIO(response.content))

plt.imshow(img)
plt.axis("off")
plt.title("Earth Image from NASA")
plt.show()
