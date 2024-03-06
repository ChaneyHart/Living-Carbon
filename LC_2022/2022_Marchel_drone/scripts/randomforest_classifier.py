import numpy as np
from sklearn.ensemble import RandomForestClassifier
from skimage import io
from osgeo import gdal

def random_forest_classification(image_path, ground_truth_path, output_path):
    # Load multiband image
    image = io.imread(image_path)

    # Load ground truth data
    ground_truth_dataset = gdal.Open(ground_truth_path)
    ground_truth = ground_truth_dataset.ReadAsArray()

    # Reshape the image and ground truth data
    rows, cols, bands = image.shape
    image_reshaped = image.reshape(rows * cols, bands)
    ground_truth_reshaped = ground_truth.reshape(rows * cols)

    # Train Random Forest classifier
    clf = RandomForestClassifier(n_estimators=100, random_state=0)
    clf.fit(image_reshaped, ground_truth_reshaped)

    # Perform classification on the image
    predicted_labels = clf.predict(image_reshaped)

    # Reshape predicted labels back to the original image shape
    predicted_labels_image = predicted_labels.reshape(rows, cols)

    # Save the classified image
    io.imsave(output_path, predicted_labels_image)
