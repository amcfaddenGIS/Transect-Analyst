"""os.environ['USE_PYGEOS'] = '0'"""
import os.path
import shutil
import sys
import io
import numpy as np
import folium
import matplotlib as plt
import rasterio
import glob
from folium import raster_layers
import pathlib as ptlib
import Transect_Analysis as TA
from PyQt5.QtWebEngineWidgets import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from PyQt5.QtWebEngineWidgets import QWebEngineView
from folium import *
import geojson
import geopandas as gpd
import folium.plugins.draw as draw
from PIL import Image


class MainWindow(QMainWindow):
    """
    The Main window that initializes the temporary directories
    The only interactive element in the main window is the resampling button
    """
    def __init__(self):
        super().__init__()
        # Initialize attributes for
        self.initUI()
        # This creates a temporary file path for rasters to be stored
        self.dir = QTemporaryDir()
        # Print temp directory for debugging
        print(self.dir.path())
        # Create path for Rasters
        raster_path = os.path.join(self.dir.path(), 'Rasters')
        # Create the raster directory and create a raster directory attribute
        os.mkdir(raster_path)
        self.raster_dir = raster_path
        # Create path for drawn trasects
        trans_path = os.path.join(self.dir.path(), 'Transects')
        # Create the drawn transect directory and create a transect directory attribute
        os.mkdir(trans_path)
        self.trans_dir = trans_path
        # Create path for sampled points, images and animations
        sample_path = os.path.join(self.dir.path(), 'Sampling')
        # Create the sampling directory and create a sampling directory attribute
        os.mkdir(sample_path)
        self.sample_dir = sample_path
        # Set the title of the windwo to transect analysis
        self.setWindowTitle("Transect Analysis")
        self.setFixedWidth(350)
        self.setFixedHeight(80)
    def initUI(self):
        # Create a QVBoxLayout
        self.raster_reproject_layout = QVBoxLayout()
        # Set the icon for the application window to a fire emoji
        self.setWindowIcon(QIcon('C:/Users/alexm/OneDrive/Desktop/GEOG 683/Test_Directory/FIRRRRRRREEEEEEEEE.png'))
        raster_reproject_form_layout = QFormLayout()
        # Create a button in the window
        self.button_2 = QPushButton('Run Raster Reprojection For Folium Map', self)
        # Fix the height and width of the buttons
        self.button_2.setFixedHeight(50)
        self.button_2.setFixedWidth(300)
        self.button_2.setStyleSheet("font-weight: bold; color: black; background-color: orange; border-color:black; border: 2px solid")
        # Connect button to function on_click
        # Clicking the button will create an output directory window
        # Clicked connect triggers functions to run
        self.button_2.clicked.connect(self.output_directory)
        raster_reproject_form_layout.addWidget(self.button_2)
        # Add a layout with the widgets
        self.raster_reproject_layout.addLayout(raster_reproject_form_layout)
        # Take the layout and set its alignment to the center of the window
        raster_reproject_form_layout.setFormAlignment(Qt.AlignCenter)
        # Create a new widget that will contain the layout and buttons
        # Set that as the central widget
        self.raster_widget = QWidget()
        self.raster_widget.setLayout(self.raster_reproject_layout)
        self.setCentralWidget(self.raster_widget)

        self.show()

    def reproject_raster(self, dir):
        """
        Function that resamples all rasters in the provided directory from the output directory function
        :param dir: Identified directory from the dialogue window
        :return: Resampled rasters in the reprojected raster folder
        """
        # Join a tif extension to the selected directory from the output directory method
        paths = os.path.join(dir + "/*.tif")
        # Identify all tif files paths and create a list object for resampling
        path = glob.glob(paths)
        # Reproject the rasters and output them to the temporary raster directory
        TA.reproject_raster(raster_path=path, output_location=self.raster_dir)

    def output_directory(self):
        """
        Open a dialog window for selecting a directory containing the tifs the user wants to resample

        :return: Directory with the tifs for resampling
        """
        # Find the current working directory
        startingDir = os.getcwd()
        # Use that directory to create a dialog window for the user to choose where their rasters are currently located
        destDir = QFileDialog.getExistingDirectory(None,
                                                   'Open Working Directory',
                                                   startingDir,
                                                   QFileDialog.ShowDirsOnly)
        if destDir is None:
            self.close()
        # The directory is then taken and the rasters with in the directory are projected to the folium CRS
        else:
            # Reproject the rasters
            self.reproject_raster(dir=destDir)
            # The open map function will create a window with a folium map that contains the reprojected raster
            self.open_map(dir=self.raster_dir)
            self.close()

    def open_map(self, dir):
        """
        Opens the folium window for transect sampling

        :param dir: Directory containing the resampled rasters
        :return: Folium Window class
        """
        # Open the folium map window
        self.w = FoliumWindow()
        # Set the folium map window directory to the reprojected raster directory
        self.w.dir = dir
        # Set the folium map windows transect directory to the transect directory from the main window
        # This is so that the drawn transect is saved to the appropriate directory
        self.w.trans_dir = self.trans_dir
        # Set the folium map windows sample directory to the sample directory from the main window
        # This is so that animations, images and sampled centroids are saved to the appropriate directory
        self.w.sample_dir = self.sample_dir
        # Show the window
        self.show()
        # Run the open map function to open the folium map
        self.w.open_map()

class FoliumWindow(QWidget):
    """
    Folium Window. Contains the folium map for transect drawing. Resampled rasters are also visualized
    for each image pass

    """
    def __init__(self):
        super().__init__()
        # All of these attributes are set before the window is shown
        self.dir = None
        self.trans_dir = None
        self.sample_dir = None
        self.raster_dir = None
        # UI attribute function
        self.initUI()
        # Title for the created window
        self.setWindowTitle("Transect Analysis")
        self.setFixedSize(1500, 800)
        self.setWindowIcon(QIcon('C:/Users/alexm/OneDrive/Desktop/GEOG 683/Test_Directory/FIRRRRRRREEEEEEEEE.png'))
        self.setStyleSheet('background-color: orange;')
        # Create a download request function to download exported transects from the Folium map.
        # This is meant to handle any download requests from websites provided in the QTWindow
        self.webview.page().profile().downloadRequested.connect(
            self.on_downloadRequested
        )
        self.show()
        # Attribute for determining the number of transects the user has drawn
        self.transect_count = 0
    def initUI(self):
        # Defined Layout
        layout = QVBoxLayout()
        map_title = QLabel()
        map_title.setText('Transect Drawing Map')
        map_title.setAlignment(Qt.AlignCenter)
        map_title.setStyleSheet("font-weight: bold; font-size: 25px; background-color: orange; border-color:black; border: 2px solid")
        layout.addWidget(map_title, 1)
        # Web view attribute that will render the folium map information
        # Interprets HTML elements
        self.webview = QWebEngineView()
        # Sets margins of the HTML map created
        # Run the open map function which adds the HTML code to the web engine view function
        # Add the web engine view to the layout
        layout.addWidget(self.webview, 9)
        box_layout = QHBoxLayout()
        # Sample Transect Button. Connects to the Sample Transects method, and runs the Sample Transects method
        # When the button is clicked
        self.sample_transects = QPushButton('Run Sampling Transect Function', self)
        self.sample_transects.setStyleSheet("font-weight: bold; color: black; background-color: orange; border-color:black; border: 2px solid")
        self.sample_transects.clicked.connect(self.Sample_Transects)
        # Transect Animation Button. Connects to the Render Animation method, and runs the Render Animation method
        # When the button is clicked.
        self.render_animation = QPushButton('Create Sampling Transect Animation', self)
        self.render_animation.setStyleSheet("font-weight: bold; color: black; background-color: orange; border-color:black; border: 2px solid")
        self.render_animation.clicked.connect(self.Render_Animation)
        # Output Animation and Images Button. Connects to the Save Animation and Images method, and runs the Save Animationa and Images method
        # When the button is clicked.
        self.export_animation = QPushButton('Output Rendered Animation to Select File Location', self)
        self.export_animation.setStyleSheet("font-weight: bold; color: black; background-color: orange; border-color:black; border: 2px solid")
        self.export_animation.clicked.connect(self.Save_Animation_and_Images)
        # Add each of the buttons to the layout
        box_layout.addWidget(self.sample_transects)
        box_layout.addWidget(self.render_animation)
        box_layout.addWidget(self.export_animation)
        # Determine how much of the window is occupied by the buttons
        layout.addLayout(box_layout, 1)
        # Set the layout
        self.setLayout(layout)
    def open_map(self):
        """
        Open the folium map in the folium window. Map will visualize the reprojected rasters, and provide the
        user the Draw plugin for drawing transects.

        :return: Folium Map
        """
        # List rasters function will pull all rasters in the temporary rasters directory
        rasters = TA.List_Rasters(self.dir)
        # Open the first raster
        r = rasterio.open(rasters[1])
        # Use the first raster to provide visualization bounds for the folium map
        min_lon, min_lat, max_lon, max_lat = r.bounds
        # Create a folium map object to visualize the reprojected rasters.
        # Bounds are determined based on the bounds of the first raster. Rather than choosing a corner of the raster, the function
        # approximates the centroid of the raster.
        m = Map(
            location=[(min_lat + max_lat)/2, (min_lon + max_lon)/2], zoom_start=13,
        )
        # For each of the rasters listed in the resampled directory, add each raster to the folium map
        for i in range(0, len(rasters)):
            # Call the first raster
            rast = rasters[i]
            # Determine the base file name of the file path for the raster
            file_name = os.path.basename(rast)
            # Get the Pass number from the base name
            name = file_name[10:17]
            # Open the raster using rasterio
            r = rasterio.open(rast)
            # Call the metadata from rasterio
            out_meta = r.meta
            # Open the first band as an array in rasterio
            r_array = r.read(1)
            # Noramlize the data so it can be easily visualized in folium.
            # If it isn't normalized the stretch becomes messed up.
            normed_data = (r_array - r_array.min()) / (r_array.max() - r_array.min())
            # Create a cmap for the reprojected raster that will be plotted
            cm = plt.cm.get_cmap('gnuplot2')
            # Use the cmap to change the color ramp of the array
            color_array = cm(normed_data)
            # Call the bounds of the raster
            min_lon, min_lat, max_lon, max_lat = r.bounds
            # Create a list of coordinates for the bounds of the raster that will be visualized in folium
            bounds = [[min_lat, min_lon], [max_lat, max_lon]]
            # Add the raster to the folium map using raster layers and imagery overlay
            raster_layers.ImageOverlay(color_array, opacity=1, bounds=bounds, name = name, show=False).add_to(m)
        # Add the ESRI satellite imagery tile to the map. Good for visualizing vegetation densities within the burn area
        # More helpful than open street map
        folium.TileLayer(
            tiles='https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}',
            attr='Esri',
            name='Esri Satellite',
            overlay=True,
            control=True
        ).add_to(m)
        # Create an HTML object that contains all HTML elements and information for the folium map
        mapObjectInHTML = m.get_name()
        # Edit the HTML object by adding javascript for outputting drawn objects by the user.
        # The javascript also includes an option for users to specify the name of the transect feature they drew, as well as a numeric ID
        m.get_root().html.add_child(folium.Element(
            # Add custom html to or javascript here #
            """
            <script type="text/javascript">
              $(document).ready(function(){           
                var x = 0;

                {map}.on("draw:created", function(e){    
                    var layer = e.layer;
                        feature = layer.feature = layer.feature || {}; 

                    var title = prompt("Please provide the name", "default");
                    var value = prompt("Please provide the value", "undefined");
                    var id = x++;

                    feature.type = feature.type || "Feature";
                    var props = feature.properties = feature.properties || {};
                    props.Id = id;
                    props.Title = title;
                    props.Value = value;
                    drawnItems.addLayer(layer);
                  });    
                });    
            </script>
            """.replace('{map}', mapObjectInHTML)))
        # Add the folium Draw plugin to the map
        folium.plugins.Draw(
            export=True,
            # This determines the position of the draw plugin
            position='topleft',
            draw_options={'polyline': {'allowIntersection': False,
                                       'transform': True},
                          'marker': {'allowIntersection': False}}
        ).add_to(m)
        folium.LayerControl().add_to(m)
        # Create a memory object to store the HTML file from the folium map object
        data = io.BytesIO()
        # Save the HTML File
        m.save(data, close_file=False)
        # Set the HTML for the webview object created earlier
        # The webviewer interprets a html file and adds the result as a widget
        self.webview.setHtml(data.getvalue().decode())
    def on_downloadRequested(self, download):
        """
        Download the drawn transect as a geojson. The download is a string containing geojson information.
        That string is converted into a geopandas dataframe and output as a geojson in the Transect temporary directory

        :param download: Downlaoded object from the download requested method
        :return: Geojson of the transect
        """
        # Read the json and convert it into a text file that will be saved into the temp directory
        geoj = download.url().path()
        # Print json to show that exists
        print(geoj)
        # Remove the excess encoding text for just the geosjon information
        Transect_Json = geoj[24:]
        # Increase the transect count attribute by 1
        self.transect_count += 1
        # Create an output location in the temporary directory for the newly created transect
        output_location = os.path.join(self.trans_dir, "Transect_{}.geojson".format(self.transect_count))
        output_location = output_location.replace(os.sep, '/')
        # Create geodataframe with the geojson string
        gdf = gpd.read_file(Transect_Json, driver = 'GeoJSON')
        # With the geodataframe, output the geojson the temporary transect folder
        gdf.to_file(output_location, driver = 'GeoJSON')
        # Create message box that lets the user know that the geojson has successfully been downloaded and transect sampling is ready
        msg = QMessageBox()
        msg.setStyleSheet("font-weight: bold; color: black; background-color: orange;")
        msg.setWindowIcon(QIcon('C:/Users/alexm/OneDrive/Desktop/GEOG 683/Test_Directory/FIRRRRRRREEEEEEEEE.png'))
        msg.setWindowTitle("Download Status")
        msg.setText("Transect {} has been downloaded and ready for sampling!".format(self.transect_count))
        msg.exec()
    def Sample_Transects(self):
        """
        Initializes a sample transect class for each raster in a user identified directory. Identified directory
        should contain the rasters the user wants to sample values from. Results of the sample transect class
        are used to create a plot showing the variation in temperature along the transect.
        Image series can be converted into an animation.


        :return: Series of Transect images for each raster in the identified directory
        """
        # Determine starting directory for dialog box
        startingDir = os.getcwd()
        # Open dialog box for the user to identify the
        rasters = QFileDialog.getExistingDirectory(None,
                                                   'Open Working Directory',
                                                   startingDir,
                                                   QFileDialog.ShowDirsOnly)
        # Call the transect geojson path
        transect = "{}/{}".format(self.trans_dir, "Transect_{}.geojson".format(self.transect_count))
        # After the transect path has been identified, use the transect path as the transect for the transect animation function
        # Outputs a series of images using the identified transect and a series of raster images
        # Full explanation of the tool is provided in the Transect Analysis Module
        TA.Transect_Animation(output_location = self.sample_dir,
                           line = transect,
                           crs=26911,
                           raster_list = rasters)
        # Provide message box to user indicating that the function has run successfully.
        msg = QMessageBox()
        msg.setStyleSheet("font-weight: bold; color: black; background-color: orange;")
        msg.setWindowIcon(QIcon('C:/Users/alexm/OneDrive/Desktop/GEOG 683/Test_Directory/FIRRRRRRREEEEEEEEE.png'))
        msg.setWindowTitle("Sample Status")
        msg.setText("Transect {} has been used for sampling!".format(self.transect_count))
        msg.exec()
        # After the user has closed the message box, the Open Results Window method is called
        self.Open_Results_Window(self.sample_dir)
    def Open_Results_Window(self, dir):
        """
        Opens a combo box window that allows users to view each output image from the sampling transect method.

        :param dir: Sampling directory containing each of the output images from the Sample Transect Method
        :return: Combo box window that contains each name of the output images
        """
        # Open the Output Plots window and use the
        self.o = Output_Plots(dir = dir)
        # Show the combo box window
        self.o.show()
    def Render_Animation(self):
        """
        Creates an animation of each output image from the sampling transect method


        :return: Gif of the sampling transect image series
        """
        # Create a directory list contains every image in the sampling directory
        self.image_dir = glob.glob("{}/*.png".format(self.sample_dir))
        # For each path in the directory list, open the image using PIL
        images = [Image.open(i) for i in self.image_dir]
        # Call the first image in the image list to be the first frame
        first_frame = images[0]
        # Create a gif by adding images to the first frame
        first_frame.save("{}/Transect_Animation.gif".format(self.sample_dir), format="GIF", append_images=images,
                       save_all=True, duration=425, loop=0)
        # Create an output message to indicate to the user that the function has finished running
        msg = QMessageBox()
        msg.setStyleSheet("font-weight: bold; color: black; background-color: orange;")
        msg.setWindowIcon(QIcon('C:/Users/alexm/OneDrive/Desktop/GEOG 683/Test_Directory/FIRRRRRRREEEEEEEEE.png'))
        msg.setWindowTitle("Animation Status")
        msg.setText("Transect Animation has been created!")
        msg.exec()
    def Save_Animation_and_Images(self):
        """
        Save each of the images and the gif animation to a directory identified by the user. Centroid and Transect shapefiles are saved as well


        :return: Saved images and gif
        """
        # Determine all file paths with the png extension sample directory
        self.image_dir = glob.glob("{}/*.png".format(self.sample_dir))
        # Determine all file paths with the gif extension in the sample directory
        gif_dir = glob.glob("{}/*.gif".format(self.sample_dir))
        # Determine the starting directory for the directory dialog window
        startingDir = os.getcwd()
        # Determine directory for outputting images and animations
        Output_Location = QFileDialog.getExistingDirectory(None,
                                                   'Open Working Directory',
                                                   startingDir,
                                                   QFileDialog.ShowDirsOnly)
        for p in self.image_dir:
            # For each file path in the image list, copy the image to the identified directory from the dialog window.
            shutil.copy(src=p, dst=Output_Location)
        # Do the same thing for the gift
        shutil.copy(gif_dir[0], dst=Output_Location)
        # Open centroids shapefile
        cent_gpd = gpd.read_file("{}/Centroids.shp".format(self.sample_dir))
        transect_gpd = gpd.read_file("{}/Transect {}.shp".format(self.sample_dir, self.transect_count))
        cent_gpd.to_file("{}/Centroids.shp".format(Output_Location))
        transect_gpd.to_file("{}/Transects.shp".format(Output_Location))
        msg = QMessageBox()
        msg.setStyleSheet("font-weight: bold; color: black; background-color: orange;")
        msg.setWindowIcon(QIcon('C:/Users/alexm/OneDrive/Desktop/GEOG 683/Test_Directory/FIRRRRRRREEEEEEEEE.png'))
        msg.setWindowTitle("Save Status")
        msg.setText("Transect Animation and Images have been saved")
        msg.exec()


class Output_Plots(QWidget):
    """
    Small window that contains a combo box. Combo box elements are the names of each output image from the sampling transect method.
    """
    def __init__(self, dir):
        super().__init__()
        self.show()
        # Flag for render image method
        self.flag = False
        self.setWindowTitle("Transect Sample Results")
        # Directory containing all image paths
        self.dir = dir
        self.setStyleSheet('background-color: orange;')
        # List containing all image paths
        self.image_dir = glob.glob("{}/*.png".format(self.dir))
        # Dictionary created by the dictionary creation method (further explanation there)
        self.path_dictionary = self.dictionary_creation()
        self.setFixedWidth(350)
        self.setFixedHeight(80)
        self.__initUI__()
    def __initUI__(self):
        # Create layout for widgets
        self.HLayout = QHBoxLayout()
        self.setWindowIcon(QIcon('C:/Users/alexm/OneDrive/Desktop/GEOG 683/Test_Directory/FIRRRRRRREEEEEEEEE.png'))
        # Create variable for combo box elements. Elements are from the Image list method (further explanation there)
        combo_box_elements = self.Image_List()
        # Combo box containing a list of image names
        self.combobox = QComboBox()
        self.combobox.setEditable(True)
        self.combobox.addItems(combo_box_elements)
        # Add reactive functionality. When the value in the combo box changes, run the render image method
        self.combobox.activated.connect(self.render_image)
        # Add the combo box widget to the layaut
        self.HLayout.addWidget(self.combobox)
        self.setLayout(self.HLayout)
    def Image_List(self):
        """
        Creates a list of image names from the image directory attribute. Each image name is added to a list and used as the
        elements for the combo box

        :return: List of image names
        """
        image_list = []
        for r in self.image_dir:
            # Get the basename from the image file path and add it to the image list
            basename = os.path.basename(r)
            image_list.append(basename)
        return image_list
    def render_image(self):
        """
        Renders an image when the user changes the combo box value. Opens a new window (render plot window)


        :return: Window with pixmap of output plot
        """

        # Call the file path value from the image name key in the path dictionary
        # i.e. of image name Pass29.png, what it's associated file path?
        # The file path is used to render a pixmap in a new plot

        image_path = self.path_dictionary[self.combobox.currentText()]
        # If the flag is still false, a new image will open and the flag will be set to true
        # If the flag is true, the current Render Plot window will close and a new window will be created
        if self.flag is False:
            # Using the image path from the path dictionary, call a window containing the image.
            self.p = Render_Plot(image = image_path)
            self.p.show()
            self.flag = True
        else:
            self.p.close()
            self.p = Render_Plot(image_path)
            self.p.show()
    def dictionary_creation(self):
        """
        Creates a dictionary containing each image name in the combo box and it's associated path.
        The dictionary is required for properly visualizing the output images from the sampling transect method


        :return: Dictionary containing file paths as values and image names as keys
        """
        # Create function that will connect the combo box text to a specific file path
        path_dict = {}
        for i in range(0, len(self.image_dir)):
            # Determine the base name of the image in the image file directory list
            base_name = os.path.basename(self.image_dir[i])
            # Create a dictionary containing the file path as the value and the image name as the key
            # This way, the image name can be used to select it's associated image path for visualization
            path_dict[base_name] = self.image_dir[i]
        return path_dict

class Render_Plot(QWidget):
    """
    Window that contains the selected image from the combo box in the Output plots window.
    """
    def __init__(self, image):
        super().__init__()
        self.setWindowTitle("Transect Sample Results")
        self.setWindowIcon(QIcon('C:/Users/alexm/OneDrive/Desktop/GEOG 683/Test_Directory/FIRRRRRRREEEEEEEEE.png'))
        self.image_path = image
        self.setFixedSize(1300, 900)
        self.show()
        self.__initUI__()
    def __initUI__(self):
        label = QLabel(self)
        # Create pix map based on the image path selected from the image dictionary
        # Image path is determined by the combo box value
        self.image = QPixmap(self.image_path)
        label.setPixmap(self.image)
        self.HLayout = QHBoxLayout()
        self.HLayout.addWidget(label)
        self.setLayout(self.HLayout)


app = QApplication(sys.argv)

window = MainWindow()
window.show()
app.exec()
window.dir.remove()
