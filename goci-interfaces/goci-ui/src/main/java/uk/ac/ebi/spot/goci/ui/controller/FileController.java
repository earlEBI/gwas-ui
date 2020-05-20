package uk.ac.ebi.spot.goci.ui.controller;

import org.apache.tomcat.util.http.fileupload.IOUtils;
import org.springframework.beans.factory.annotation.Value;
import org.springframework.core.io.Resource;
import org.springframework.http.HttpStatus;
import org.springframework.http.MediaType;
import org.springframework.stereotype.Controller;
import org.springframework.web.bind.annotation.ExceptionHandler;
import org.springframework.web.bind.annotation.RequestMapping;
import org.springframework.web.bind.annotation.RequestMethod;
import org.springframework.web.bind.annotation.ResponseBody;
import org.springframework.web.bind.annotation.ResponseStatus;

import javax.servlet.http.HttpServletResponse;
import java.io.*;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Properties;

/**
 * Created by emma on 24/02/15.
 *
 * @author emma
 *         <p>
 *         Controller used to handle download of files and generation of stats
 */
@Controller
public class FileController {

    // These parameters are read from application.properties file
    @Value("${download.full}")
    private Resource fullFileDownload;

    @Value("${download.alternative}")
    private Resource alternativeFileDownload;

    @Value("${download.studies}")
    private Resource studiesFileDownload;

    @Value("${download.studiesAlternative}")
    private Resource alternativeStudiesDownload;

    @Value("${download.ancestry}")
    private Resource ancestryFileDownload;

    @Value("${download.efoMappings}")
    private Resource efoMappingsDownload;

    @Value("${download.NCBI}")
    private Resource fullFileDownloadNcbi;

    @Value("${catalog.stats.file}")
    private Resource catalogStatsFile;

    @Value("${summary.stats.file}")
    private Resource summaryStatsFile;

//    @Value("${summary.stats.fullpvalue.file}")
//    private Resource summaryStatsFullPValueFile;
//
    @Value("${download.unpublished.studies}")
    private Resource unpublishedStudiesFileDownload;

    @Value("${download.unpublished.ancestries}")
    private Resource unpublishedAncestriesFileDownload;

//    @Value("${download.ensemblmapping}")
//    private Resource ensemblMappingFileDownload;

    @RequestMapping(value = "api/search/downloads/full",
                    method = RequestMethod.GET)
    public void getFullDownload(HttpServletResponse response) throws IOException {
        if (fullFileDownload.exists() && catalogStatsFile.exists()) {

            Properties properties = new Properties();
            properties.load(catalogStatsFile.getInputStream());
            String releasedate = properties.getProperty("releasedate");
            String ensemblbuild = properties.getProperty("ensemblbuild");


            String fileName = "gwas_catalog_v1.0-associations_e".concat(ensemblbuild).concat("_r").concat(releasedate).concat(".tsv");
            buildDownload(fileName, fullFileDownload.getInputStream(), response);
        }
        else {
            throw new FileNotFoundException();
        }
    }

//    @RequestMapping(value = "api/downloads/fullpvalue",
//            method = RequestMethod.GET, produces = MediaType.APPLICATION_JSON_VALUE)
//    public void getFullPValueStudiesDownload(HttpServletResponse response) throws IOException {
//        String responseString = null;
//        if (summaryStatsFullPValueFile.exists()) {
//            byte[] bytes = Files.readAllBytes(summaryStatsFullPValueFile.getFile().toPath());
//            IOUtils.copy(new BufferedInputStream(new ByteArrayInputStream(bytes)),
//                    new BufferedOutputStream(response.getOutputStream()));
////            buildJsonDownload(summaryStatsFullPValueFile.getInputStream(), response);
//            responseString = new String(bytes);
//        }
//        else {
//            throw new FileNotFoundException();
//        }
//    }

    @RequestMapping(value = "api/search/downloads/unpublished_studies",
            method = RequestMethod.GET)
    public void getUnpublishedStudiesDownload(HttpServletResponse response) throws IOException {
        if (unpublishedStudiesFileDownload.exists() && catalogStatsFile.exists()) {

            Properties properties = new Properties();
            properties.load(catalogStatsFile.getInputStream());
            String releasedate = properties.getProperty("releasedate");

            String fileName = "gwas-catalog-unpublished-studies-r".concat(releasedate).concat("-v1.0.3.tsv");
            buildDownload(fileName, unpublishedStudiesFileDownload.getInputStream(), response);
        }
        else {
            throw new FileNotFoundException();
        }
    }

    @RequestMapping(value = "api/search/downloads/unpublished_ancestries",
            method = RequestMethod.GET)
    public void getUnpublishedAncestriesDownload(HttpServletResponse response) throws IOException {
        if (unpublishedAncestriesFileDownload.exists() && catalogStatsFile.exists()) {

            Properties properties = new Properties();
            properties.load(catalogStatsFile.getInputStream());
            String releasedate = properties.getProperty("releasedate");

            String fileName = "gwas-catalog-unpublished-studies-r".concat(releasedate).concat("-v1.0.3.tsv");
            buildDownload(fileName, unpublishedAncestriesFileDownload.getInputStream(), response);
        }
        else {
            throw new FileNotFoundException();
        }
    }

    @RequestMapping(value = "api/search/downloads/studies",
                    method = RequestMethod.GET)
    public void getStudiesDownload(HttpServletResponse response) throws IOException {
        if (studiesFileDownload.exists() && catalogStatsFile.exists()) {

            Properties properties = new Properties();
            properties.load(catalogStatsFile.getInputStream());
            String releasedate = properties.getProperty("releasedate");

            String fileName = "gwas_catalog_v1.0-studies_r".concat(releasedate).concat(".tsv");
            buildDownload(fileName, studiesFileDownload.getInputStream(), response);
        }
        else {
            throw new FileNotFoundException();
        }
    }


    @RequestMapping(value = "api/search/downloads/alternative",
                    method = RequestMethod.GET,
                    produces = MediaType.TEXT_PLAIN_VALUE)
    public void getAlternativeDownload(HttpServletResponse response) throws IOException {
        if (alternativeFileDownload.exists() && catalogStatsFile.exists()) {

            Properties properties = new Properties();
            properties.load(catalogStatsFile.getInputStream());
            String releasedate = properties.getProperty("releasedate");
            String ensemblbuild = properties.getProperty("ensemblbuild");

            String fileName = "gwas_catalog_v1.0.2-associations_e".concat(ensemblbuild).concat("_r").concat(releasedate).concat(".tsv");
            buildDownload(fileName, alternativeFileDownload.getInputStream(), response);
        }
        else {
            throw new FileNotFoundException();
        }
    }

    @RequestMapping(value = "api/search/downloads/studies_alternative",
                    method = RequestMethod.GET)
    public void getAlternativeStudiesDownload(HttpServletResponse response) throws IOException {
        if (alternativeStudiesDownload.exists() && catalogStatsFile.exists()) {

            Properties properties = new Properties();
            properties.load(catalogStatsFile.getInputStream());
            String releasedate = properties.getProperty("releasedate");


            String fileName = "gwas_catalog_v1.0.2-studies_r".concat(releasedate).concat(".tsv");
            buildDownload(fileName, alternativeStudiesDownload.getInputStream(), response);
        }
        else {
            throw new FileNotFoundException();
        }
    }

    @RequestMapping(value = "api/search/downloads/trait_mappings",
                    method = RequestMethod.GET)
    public void getTraitMappingsDownload(HttpServletResponse response) throws IOException {
        if (efoMappingsDownload.exists() && catalogStatsFile.exists()) {

            Properties properties = new Properties();
            properties.load(catalogStatsFile.getInputStream());
            String releasedate = properties.getProperty("releasedate");


            String fileName = "gwas_catalog_trait-mappings_r".concat(releasedate).concat(".tsv");
            buildDownload(fileName, efoMappingsDownload.getInputStream(), response);
        }
        else {
            throw new FileNotFoundException();
        }
    }


    @RequestMapping(value = "api/search/downloads/full_NCBI",
                    method = RequestMethod.GET,
                    produces = MediaType.TEXT_PLAIN_VALUE)
    public void getFullNcbiDownload(HttpServletResponse response) throws IOException {

        if (fullFileDownloadNcbi.exists()) {

            InputStream inputStream = null;
            inputStream = fullFileDownload.getInputStream();

            OutputStream outputStream;
            outputStream = response.getOutputStream();

            IOUtils.copy(new BufferedInputStream(inputStream), new BufferedOutputStream(outputStream));
            inputStream.close();
            outputStream.close();

        }
        else {
            throw new FileNotFoundException();
        }

    }

    @RequestMapping(value = "api/search/stats", method = RequestMethod.GET, produces = "application/json")
    public @ResponseBody Map<String, Object> getCatalogStats() {
        Map<String, Object> response = new HashMap<>();

        String releasedate;
        String studycount;
        String snpcount;
        String associationcount;
        String genebuild;
        String dbsnpbuild;
        String ensemblbuild;

        Properties properties = new Properties();
        try {
            properties.load(catalogStatsFile.getInputStream());
            releasedate = properties.getProperty("releasedate");
            studycount = properties.getProperty("studycount");
            snpcount = properties.getProperty("snpcount");
            associationcount = properties.getProperty("associationcount");
            genebuild = properties.getProperty("genomebuild");
            dbsnpbuild = properties.getProperty("dbsnpbuild");
            ensemblbuild = properties.getProperty("ensemblbuild");

            response.put("date", releasedate);
            response.put("studies", studycount);
            response.put("snps", snpcount);
            response.put("associations", associationcount);
            response.put("genebuild", genebuild);
            response.put("dbsnpbuild", dbsnpbuild);
            response.put("ensemblbuild", ensemblbuild);

        }
        catch (IOException e) {
            throw new RuntimeException(
                    "Unable to return catolog stats: failed to read catalog.stats.file resource", e);
        }

        return response;
    }

    @RequestMapping(value = "api/search/summaryStatsResources", method = RequestMethod.GET, produces = "application/json")
    public @ResponseBody Map<String, Object> getSummaryStatsResources() {
        Map<String, Object> response = new HashMap<>();

        List<String> resources = new ArrayList<>();

        try {
            if(summaryStatsFile.exists()){
                InputStream in = new BufferedInputStream(summaryStatsFile.getInputStream());
                BufferedReader reader = new BufferedReader(new InputStreamReader(in));
                String line;
                while ((line = reader.readLine()) != null) {
                    resources.add(line);
                }
                in.close();
                reader.close();
                response.put("resources", resources);
            }
        }
        catch (IOException e) {
            throw new RuntimeException(
                    "Unable to return summary stats resources: failed to read summary.stats.file resource", e);
        }

        return response;
    }


    @RequestMapping(value = "api/search/downloads/ancestry",
                    method = RequestMethod.GET)
    public void getAncestryDownload(HttpServletResponse response) throws IOException {
        if (ancestryFileDownload.exists() && catalogStatsFile.exists()) {

            Properties properties = new Properties();
            properties.load(catalogStatsFile.getInputStream());
            String releasedate = properties.getProperty("releasedate");

            String fileName = "gwas_catalog-ancestry_r".concat(releasedate).concat(".tsv");
            buildDownload(fileName, ancestryFileDownload.getInputStream(), response);

        }
        else {
            throw new FileNotFoundException();
        }
    }


    @ResponseStatus(value = HttpStatus.NOT_FOUND, reason = "File not found for download")
    @ExceptionHandler(FileNotFoundException.class)
    public void FileNotFoundException(FileNotFoundException fileNotFoundException) {
    }

    private void buildDownload(String fileName, InputStream inputStream, HttpServletResponse response)
            throws IOException {
        response.setContentType("text/tsv");
        response.setHeader("Content-Disposition", "attachement; filename=" + fileName);

        OutputStream outputStream;
        outputStream = response.getOutputStream();

        IOUtils.copy(new BufferedInputStream(inputStream), new BufferedOutputStream(outputStream));
        inputStream.close();
        outputStream.close();

    }
}
