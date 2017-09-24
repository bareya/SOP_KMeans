# SOP K-Means

Recently I faced a need of slightly more efficient approach to compute clusters from a point cloud.

Cluster SOP seems to be a great tool, but unfortunately is a bit slow, and it doesn't take the advantage of all cpu cores.

Here you can find an attemp to make cluster calculation a bit faster. Code was used for different purpouses, but to make your life more easier I wrapped it as a SOP operator.  

![stats](https://github.com/bareya/SOP_KMeans/blob/master/images/stats.png)

`k-means++` is not implemented yet, I might consider implementing it in the future.
