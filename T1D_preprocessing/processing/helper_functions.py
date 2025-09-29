import numpy as np
import cv2 as cv
import scipy.optimize
import sklearn.neighbors
import time
import sys
from typing import Optional
import matplotlib.pyplot as plt
import seaborn as sns

def icp(a, b,
        max_iter = 10,
        max_time:Optional[int] = 4,
        ):
    """
    Iterative Closest Point (ICP) algorithm
    
    Parameters
    ----------
    a
        dataset a, whose shape is [n,2]
    b
        dataset a, whose shape is [m,2]
    max_time
        max iter times
    
    Return
    ----------
    T_opt
        optimal transform matrix, please see https://docs.opencv.org/3.4/da/d6e/tutorial_py_geometric_transformations.html
    error_max
        error
    
    Refer
    ----------
    (https://github.com/gao-lab/SLAT)
    (https://github.com/gao-lab/SLAT/tree/main/scSLAT/model/prematch)
    https://stackoverflow.com/questions/20120384/iterative-closest-point-icp-implementation-on-python
    """
    # Calculate residual errors between src and dst points after applying transformation
    def res(p,src,dst):
        T = np.matrix([[np.cos(p[2]),-np.sin(p[2]),p[0]],
        [np.sin(p[2]), np.cos(p[2]),p[1]],
        [0 ,0 ,1 ]])
        n = np.size(src,0)
        xt = np.ones([n,3])
        xt[:,:-1] = src
        xt = (xt*T.T).A
        d = np.zeros(np.shape(src))
        d[:,0] = xt[:,0]-dst[:,0]
        d[:,1] = xt[:,1]-dst[:,1]
        r = np.sum(np.square(d[:,0])+np.square(d[:,1]))
        return r

    # Calculate Jacobian matrix
    def jac(p,src,dst):
        T = np.matrix([[np.cos(p[2]),-np.sin(p[2]),p[0]],
        [np.sin(p[2]), np.cos(p[2]),p[1]],
        [0 ,0 ,1 ]])
        n = np.size(src,0)
        xt = np.ones([n,3])
        xt[:,:-1] = src
        xt = (xt*T.T).A
        d = np.zeros(np.shape(src))
        d[:,0] = xt[:,0]-dst[:,0]
        d[:,1] = xt[:,1]-dst[:,1]
        dUdth_R = np.matrix([[-np.sin(p[2]),-np.cos(p[2])],
                            [ np.cos(p[2]),-np.sin(p[2])]])
        dUdth = (src*dUdth_R.T).A
        g = np.array([  np.sum(2*d[:,0]),
                        np.sum(2*d[:,1]),
                        np.sum(2*(d[:,0]*dUdth[:,0]+d[:,1]*dUdth[:,1])) ])
        return g
    
    def hess(p,src,dst):
        n = np.size(src,0)
        T = np.matrix([[np.cos(p[2]),-np.sin(p[2]),p[0]],
        [np.sin(p[2]), np.cos(p[2]),p[1]],
        [0 ,0 ,1 ]])
        n = np.size(src,0)
        xt = np.ones([n,3])
        xt[:,:-1] = src
        xt = (xt*T.T).A
        d = np.zeros(np.shape(src))
        d[:,0] = xt[:,0]-dst[:,0]
        d[:,1] = xt[:,1]-dst[:,1]
        dUdth_R = np.matrix([[-np.sin(p[2]),-np.cos(p[2])],[np.cos(p[2]),-np.sin(p[2])]])
        dUdth = (src*dUdth_R.T).A
        H = np.zeros([3,3])
        H[0,0] = n*2
        H[0,2] = np.sum(2*dUdth[:,0])
        H[1,1] = n*2
        H[1,2] = np.sum(2*dUdth[:,1])
        H[2,0] = H[0,2]
        H[2,1] = H[1,2]
        d2Ud2th_R = np.matrix([[-np.cos(p[2]), np.sin(p[2])],[-np.sin(p[2]),-np.cos(p[2])]])
        d2Ud2th = (src*d2Ud2th_R.T).A
        H[2,2] = np.sum(2*(np.square(dUdth[:,0])+np.square(dUdth[:,1]) + d[:,0]*d2Ud2th[:,0]+d[:,0]*d2Ud2th[:,0]))
        return H
        
    
    assert a.shape[0] == 2 and b.shape[0] == 2, "a and b must be 2xn"
    t0 = time.time()
    init_pose = (0,0,0)
    src = np.array([a.T], copy=True).astype(np.float32)
    dst = np.array([b.T], copy=True).astype(np.float32)
    Tr = np.array([[np.cos(init_pose[2]),-np.sin(init_pose[2]),init_pose[0]],
                   [np.sin(init_pose[2]), np.cos(init_pose[2]),init_pose[1]],
                   [0,                    0,                   1          ]]).astype(np.float32)
    #print("src",np.shape(src))
    #print("Tr[0:2]",np.shape(Tr[0:2]))
    src = cv.transform(src, Tr[0:2,:])
    p_opt = np.array(init_pose)
    T_opt = np.array([])
    error_max = sys.maxsize
    first = False
    iteration = 0
    indices = []
    #while not(first and time.time() - t0 > max_time):
    while (iteration < max_iter and (time.time() - t0 < max_time)):
        iteration += 1
        distances, indices = sklearn.neighbors.NearestNeighbors(n_neighbors=1, algorithm='auto',p = 3).fit(dst[0]).kneighbors(src[0])
        p = scipy.optimize.minimize(res,[0,0,0],args=(src[0],dst[0, indices.T][0]),method='Newton-CG',jac=jac,hess=hess).x
        T  = np.array([[np.cos(p[2]),-np.sin(p[2]),p[0]],[np.sin(p[2]), np.cos(p[2]),p[1]]])
        p_opt[:2]  = (p_opt[:2]*np.matrix(T[:2,:2]).T).A       
        p_opt[0] += p[0]
        p_opt[1] += p[1]
        p_opt[2] += p[2]
        #src = cv.warpAffine(src, T, (src.shape[1], src.shape[0]))
        src = cv.transform(src, Tr[0:2,:])
        Tr = (np.matrix(np.vstack((T,[0,0,1])))*np.matrix(Tr)).A
        error = res([0,0,0],src[0],dst[0, indices.T][0])

        if error < error_max:
            error_max = error
            first = True
            T_opt = Tr
            indices_opt = indices

    p_opt[2] = p_opt[2] % (2*np.pi)
    
    # return optimal transform matrix, the max error, and the indices in b that correspond to a
    return T_opt, error_max, indices_opt


## Pad images to the same size
def pad_images(img1_padded, image2, pad_color):

    # Image shapes
    y1, x1 = np.squeeze(img1_padded).shape
    y2, x2 = np.squeeze(image2).shape
    
    max_width = np.maximum(x1, x2)
    max_height = np.maximum(y1, y2)

    # Compute offsets between image centers
    offset_x1 = max((x2 - x1) // 2, 0)
    offset_y1 = max((y2 - y1) // 2, 0)
    offset_x2 = max((x1 - x2) // 2, 0)
    offset_y2 = max((y1 - y2) // 2, 0)

    # Create empty padded images
    img1_padded_padded = np.full((max_height, max_width), pad_color, dtype=np.uint8)
    image2_padded = np.full((max_height, max_width), pad_color, dtype=np.uint8)

    # Copy original images into the center of the padded images
    img1_padded_padded[offset_y1 : offset_y1 + y1,
                  offset_x1 : offset_x1 + x1] = img1_padded.copy()
                                                                         
    image2_padded[offset_y2 : offset_y2 + y2,
                  offset_x2 : offset_x2 + x2] = image2.copy()
    
    # Make tupple with offset coordinates to return
    img1_offset_coord = (offset_y1, offset_y1+y1, offset_x1, offset_x1+x1)
    img2_offset_coord = (offset_y2, offset_y2+y2, offset_x2, offset_x2+x2)
    
    return img1_padded_padded, image2_padded, img1_offset_coord, img2_offset_coord

## Plot the centroids of the islets
def plot_centroids(ax, df, title=None):
    sns.scatterplot(x="centroid_1_dx", y="centroid_0_dy", s=20, hue="panel_islet_parent", data=df, ax=ax, alpha = 0.7)
    ax.set_xlabel("centroid1: X-Coordinate")
    ax.set_ylabel("Centroid0: Y-Coordinate")
    ax.legend(title=title)
    ax.invert_yaxis()

## Write out the dataframes to .csv files.
def write_out_dfs(df, panel_names, idx, df_post_icp, inter_panel_neighbors):

    df_post_icp["centroid-1"] = df_post_icp["centroid_1_dx"]
    df_post_icp["centroid-0"] = df_post_icp["centroid_0_dy"]

    ## Write out the dataframes to .csv files.
    for panel_name in panel_names:
        seg_subdir = "seg_subdir" + "_" + panel_name
        img_path = "img_path" + "_" + panel_name
        df_temp = df_post_icp[df_post_icp["Panel"] == panel_name].copy()
        df_temp = df_temp.drop(columns=["Panel", "panel_islet_parent", "centroid_1_dx", "centroid_0_dy"])
        
        ## Registered regionprops:
        isl_regionprops_reg_file = df.loc[:,(seg_subdir, "registered_regionprops")].iloc[idx] / df.loc[:,(img_path, "registered_regionprops")].iloc[idx]
        df_temp.to_csv(str(isl_regionprops_reg_file))
        print(isl_regionprops_reg_file)
        
        ## Inter-panel neighbors:
        ## Note this is the same file:
        isl_regionprops_nei_file = df.loc[:,(seg_subdir, "inter_panel_neighbors")].iloc[idx] / df.loc[:,(img_path, "inter_panel_neighbors")].iloc[idx]
        inter_panel_neighbors["Islet_cell_id"] = inter_panel_neighbors["islet_roi_id"].astype(str) + "_" + inter_panel_neighbors["Object_Islet"].astype(str)
        inter_panel_neighbors["Immune_cell_id"] = inter_panel_neighbors["immune_roi_id"].astype(str) + "_" + inter_panel_neighbors["Object_Immune"].astype(str)
        inter_panel_neighbors.to_csv(str(isl_regionprops_nei_file))
        print(isl_regionprops_nei_file)
