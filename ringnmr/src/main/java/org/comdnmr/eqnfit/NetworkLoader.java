package org.comdnmr.eqnfit;

import java.io.File;
import java.io.IOException;

import java.nio.file.Files;
import java.nio.file.StandardCopyOption;

import java.util.HashMap;
import java.util.Map;

import org.tensorflow.SavedModelBundle;

import org.nmrfx.utilities.UnZipper;

/**
 *
 * @author simonhulse
 */
class NetworkLoader {
    private static NetworkLoader instance = null;

    private Map<String, SavedModelBundle> loadedNetworks;

    private File zipSrc = new File("/data/parameter-networks.zip");
    private File zipDst = new File(System.getProperty("java.io.tmpdir"), "parameter-networks.zip");
    private File zipDstParent = zipDst.getParentFile();
    private File rootDir = new File(zipDstParent, "parameter-networks");

    private NetworkLoader() throws IOException {
        if (!zipDst.exists()) {
            // TODO: Gracfully handle IOException
            Files.copy(
                getClass().getResourceAsStream(zipSrc.getPath()),
                zipDst.toPath(),
                StandardCopyOption.REPLACE_EXISTING);
            UnZipper unZipper = new UnZipper(zipDstParent, zipDst.getPath());
            unZipper.unzip();
        }

        loadedNetworks = new HashMap<>();
    }

    public static synchronized NetworkLoader getNetworkLoader() throws IOException {
        if (instance == null) {
            instance = new NetworkLoader();
        }

        return instance;
    }

    public SavedModelBundle fetchNetwork(String path) {
        SavedModelBundle network;
        if (loadedNetworks.containsKey(path)) {
            network = loadedNetworks.get(path);
        } else {
            File absolutePath = getAbsolutePath(path);
            network = SavedModelBundle.load(absolutePath.getPath(), "serve");
            loadedNetworks.put(path, network);
        }
        return network;
    }

    private File getAbsolutePath(String relativePath) {
        return new File(rootDir, relativePath);
    }
}
