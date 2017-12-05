## Class Server

Server: `gpuvannberg.biology.gatech.edu`
Data Folder: `/data2/AHCG2017FALL/`

# Virtual Machine

Prebuilt VMs were stored in `/data2/VMbox_prebuilt`. 

  Monaco, Christopher                    cmonaco3         Ubuntu-64-DR-AHCG2017-p10025.ova    10025
  
## Adding and accessing VM

VM was imported with this command
 
    $ vboxmanage import Ubuntu-64-DR-AHCG2017-p10025.ova

To get list of VMs to find correct ID

    $ vboxmanage list vms

VM was started with

    $ vboxmanage startvm Ubuntu-64-DR-AHCG2017 --type headless

VM was logged into with

    $ ssh vannberglab@localhost -p 10025
    password: vanberglab
  
WM can be shutdown with

    $ vboxmanage controlvm Ubuntu-64-DR-AHCG2017 poweroff soft

## Copying files from GPUVannberg to VM

On GPUVannberg:

    $ scp -r -P 10025 /data2/AHCG2017FALL/bin/ vannebrglab@localhost:~/
  
On VM:
  
    $ scp -r cmonaco3@gpuvannberg.biology.gatech.edu:/data2/AHCG2017FALL/data/ .

## Cloning and Increasing Disk Size

Cloning disk

    $ vboxmanage clonehd Ubuntu-64-DR-AHCG2017-p10025-disk001.vmdk Ubuntu-64-DR-AHCG2017.vdi --format vdi
    
Resizing disk to 120 gig

    $  vboxmanage modifyhd Ubuntu-64-DR-AHCG2017.vdi --resize 120000


# Read Depth Coverage Calculations

## Extract Regions of Interest from BAM File

Extracting only certain regions significantly increases search time.

   $ samtools view input.bam "Chr10:18000-45500" > output.bam
   
## Compute Genome Coverage using Bed Tools

We will use the genomecov tool. [documentation](http://bedtools.readthedocs.io/en/stable/content/tools/genomecov.html)

   $ bedtools genomecov -d -ibam input.bam -g genome.bed > output.coverage
   
There's a script in bin/scripts
