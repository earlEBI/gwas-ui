package uk.ac.ebi.spot.goci.ui.controller;

import com.fasterxml.jackson.databind.JsonNode;
import com.fasterxml.jackson.databind.ObjectMapper;
import com.fasterxml.jackson.databind.node.*;
import com.fasterxml.jackson.databind.util.JSONPObject;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.http.MediaType;
import org.springframework.http.ResponseEntity;
import org.springframework.http.converter.json.GsonBuilderUtils;
import org.springframework.stereotype.Controller;
import org.springframework.ui.Model;
import org.springframework.web.bind.annotation.PathVariable;
import org.springframework.web.bind.annotation.RequestMapping;
import org.springframework.web.bind.annotation.RequestParam;
import org.springframework.web.bind.annotation.ResponseBody;
import org.springframework.web.client.RestTemplate;
import uk.ac.ebi.spot.goci.model.SnpResult;
import uk.ac.ebi.spot.goci.ui.SearchConfiguration;

import java.io.IOException;
import java.lang.reflect.Array;
import java.util.*;

/**
 * Created by laurent on 12/01/15.
 */

@Controller
public class SNPController {

    private SearchConfiguration searchConfiguration;

    @Autowired
    private ObjectMapper mapper;

    @Autowired
    public SNPController(SearchConfiguration searchConfiguration) {
        this.searchConfiguration = searchConfiguration;
    }

    @RequestMapping(value = "/variants/{rsId}", produces = MediaType.TEXT_HTML_VALUE) String search(Model model,
                                                                                               @PathVariable String rsId,
                                                                                               @RequestParam(required = false) String filter) {
        SnpResult result = new SnpResult();
        result.setQuery(rsId);
        result.setFilter(filter);
        result.setRsId(rsId);
        model.addAttribute("result", result);
        return "/variant-page";
    }

    @RequestMapping(value = "/variantData/{rsId}")
    @ResponseBody
    Map<String, Object> variantData(@PathVariable String rsId) {

        Map<String, Object> map = new HashMap<>();
        try {//http://localhost:8080/api/search?q=rsID%3A%22rs7329174%22+AND+resourcename%3Avariant
            RestTemplate template = new RestTemplate();

            String url = "http://localhost:8080/api/search?q=rsID:\"" + rsId  + "\" AND resourcename:variant";
            String response = template.getForObject(url, String.class);
            JsonNode solrNode = mapper.readTree(response);
            JsonNode node =
                    template.getForObject("http://localhost:8081/api/singleNucleotidePolymorphisms/" + rsId +
                            "/associations?projection=associationByStudy", JsonNode.class);
            Set<String> mappedGenes = new HashSet<>();
            for(JsonNode doc: solrNode.get("response").get("docs")){
                for(JsonNode geneNode: doc.get("mappedGenes")){
                    String mappedGene = geneNode.asText().split("\\|")[0];
                    mappedGenes.add(mappedGene);
                }
            }
            Map<String, JsonNode> efoMap = new TreeMap<>();
            JsonNode efoNode = buildEfoNode(rsId, node, efoMap);
            map.put("efo", efoNode);
            Map<String, JsonNode> studyMap = new TreeMap<>();
            JsonNode studyNode = buildStudyNode(rsId, node, studyMap);
            map.put("studies", studyNode);
            node =
                    template.getForObject("http://localhost:8081/api/singleNucleotidePolymorphisms/" + rsId +
                            "/associations?projection=associationBySnp", JsonNode.class);
            map.put("associations", buildAssociationNode(rsId, node, efoMap, studyMap, mappedGenes));
            node = template.getForObject("http://localhost:8081/api/singleNucleotidePolymorphisms/" + rsId,
                    JsonNode.class);
            map.put("snp", node);
        } catch (Exception e) {
            e.printStackTrace();
        }
        return map;
    }

    private JsonNode buildStudyNode(String rsId, JsonNode node, Map<String, JsonNode> studyMap) {
        ArrayNode studyNode = JsonNodeFactory.instance.arrayNode();
        for(JsonNode associationObject : node.get("_embedded").get("associations")) {
            ObjectNode innerNode = studyNode.addObject();
            String associationUrl = associationObject.get("_links").get("self").get("href").asText();
            String associationId = associationUrl.substring(associationUrl.lastIndexOf('/') + 1);
            JsonNode stNode = associationObject.get("study");
            studyMap.put(associationId, stNode);
            String accessionId = stNode.get("accessionId").asText();
            innerNode.put("id", accessionId);
            innerNode.put("accessionId", accessionId);
            ArrayNode gtNode = innerNode.putArray("genotypingTechnologies");
            for(JsonNode tech: stNode.get("genotypingTechnologies")){
                gtNode.add(tech.get("genotypingTechnology"));
            }
            innerNode.put("fullPvalueSet", stNode.get("fullPvalueSet").asText());
            innerNode.put("author_s", stNode.get("publicationInfo").get("author").get("fullname").asText());
            RestTemplate template = new RestTemplate();
            JsonNode assocNode =
                    template.getForObject("http://localhost:8081/api/studies/" + accessionId +
                            "/associations", JsonNode.class);
            //get count via /api/studies/GCST001785/associations
            ArrayNode assocArray = (ArrayNode) assocNode.get("_embedded").get("associations");
            innerNode.put("associationCount", assocArray.size());
            innerNode.put("publicationDate", stNode.get("publicationInfo").get("publicationDate").asText());
            innerNode.put("publication", stNode.get("publicationInfo").get("publication").asText());
            innerNode.put("title", stNode.get("publicationInfo").get("title").asText());
            innerNode.put("pubmedId", stNode.get("publicationInfo").get("pubmedId").asText());
            innerNode.put("traitName_s", stNode.get("diseaseTrait").get("trait").asText());
            innerNode.set("efoLink", addEfoLinks(associationObject));
            ArrayNode ancestryNode = innerNode.putArray("ancestryLinks");
            for(JsonNode ancestry: stNode.get("ancestries")){
                String type = ancestry.get("type").asText();
                String size = ancestry.get("numberOfIndividuals").asText();
                ArrayNode countryNode = (ArrayNode) ancestry.get("countryOfRecruitment");
                String country = countryNode == null ? "" : countryNode.get(0).get("countryName").asText();
                ArrayNode groupNode = (ArrayNode) ancestry.get("ancestralGroups");
                String group = groupNode == null ? "" : groupNode.get(0).get("ancestralGroup").asText();
                String link = type + "||" + country + "|" + group  + "|" + size + "|";
                ancestryNode.add(link);
            }
            innerNode.put("initialSampleDescription", stNode.get("initialSampleSize").asText());
            innerNode.put("replicateSampleDescription", stNode.get("replicationSampleSize").asText());
            //ancestrylinks "replication|NR|Thailand|South East Asian|1413|Bangkok, Thailand"

        }
        return studyNode;
    }

        private JsonNode buildEfoNode(String rsId, JsonNode node, Map<String, JsonNode> efoMap){
        ArrayNode efoNode = JsonNodeFactory.instance.arrayNode();
        for(JsonNode associationObject : node.get("_embedded").get("associations")){
            String associationUrl = associationObject.get("_links").get("self").get("href").asText();
            String associationId = associationUrl.substring(associationUrl.lastIndexOf('/') + 1);
            //get efos,
            ObjectNode asssociationNode = JsonNodeFactory.instance.objectNode();
            String diseaseTrait = associationObject.get("study").get("diseaseTrait").get("trait").textValue();
            asssociationNode.put("diseaseTrait", diseaseTrait);
            for(JsonNode innerEfoCode: associationObject.get("efoTraits")){
                String shortForm = innerEfoCode.get("shortForm").textValue();
                String traitName = innerEfoCode.get("trait").textValue();
                String uri = innerEfoCode.get("uri").textValue();
                ObjectNode innerNode = JsonNodeFactory.instance.objectNode();
                innerNode.set("shortForm", JsonNodeFactory.instance.arrayNode().add(shortForm));
                innerNode.set("mappedLabel", JsonNodeFactory.instance.arrayNode().add(traitName));
                innerNode.set("traitName", JsonNodeFactory.instance.arrayNode().add(diseaseTrait));
                innerNode.set("mappedUri", JsonNodeFactory.instance.arrayNode().add(uri));
                innerNode.set("rsId", JsonNodeFactory.instance.arrayNode().add(rsId));
                efoNode.add(innerNode);
                asssociationNode.set(rsId, innerNode);
            }
            efoMap.put(associationId, asssociationNode);
        }
        return efoNode;
    }

    private JsonNode buildAssociationNode(String rsId, JsonNode node, Map<String, JsonNode> efoMap,
                                          Map<String, JsonNode> studyMap, Set<String> mappedGenes){
        ArrayNode associationArrayNode = JsonNodeFactory.instance.arrayNode();
        for(JsonNode associationObject : node.get("_embedded").get("associations")){
            String associationUrl = associationObject.get("_links").get("self").get("href").asText();
            String associationId = associationUrl.substring(associationUrl.lastIndexOf('/') + 1);
            String diseaseTrait = efoMap.get(associationId).get("diseaseTrait").asText();
            ObjectNode associationNode = JsonNodeFactory.instance.objectNode();
            ArrayNode strongestRiskNode = JsonNodeFactory.instance.arrayNode();
            for(JsonNode lociNode : associationObject.get("loci")){
                for(JsonNode allele: lociNode.get("strongestRiskAlleles")){
                    strongestRiskNode.add(allele.get("riskAlleleName").textValue());
                }
            }
            associationNode.set("strongestAllele", strongestRiskNode);
            associationNode.put("riskFrequency", associationObject.get("riskFrequency").textValue());
            associationNode.put("pValueMantissa", associationObject.get("pvalueMantissa").asText());
            associationNode.put("pValueExponent", associationObject.get("pvalueExponent").asText());
            if(associationObject.get("pvalueDescription") != null && !associationObject.get("pvalueDescription").isNull()) {
                associationNode.putArray("qualifier").add(associationObject.get("pvalueDescription").asText());
            }

            associationNode.put("orPerCopyNum", associationObject.get("orPerCopyNum").asText());
            if(associationObject.get("orDescription") != null && !associationObject.get("orDescription").isNull()){
                associationNode.put("orDescription", associationObject.get("orDescription").textValue());
            }
            if(associationObject.get("betaNum") != null && !associationObject.get("betaNum").isNull()) {
                associationNode.put("betaNum", associationObject.get("betaNum").asText());
            }
            associationNode.put("betaUnit", associationObject.get("betaUnit").asText());
            associationNode.put("betaDirection", associationObject.get("betaDirection").asText());
            associationNode.put("range", associationObject.get("range").textValue());

            ArrayNode mappedGeneNode = JsonNodeFactory.instance.arrayNode();
            for(String geneName: mappedGenes){
                mappedGeneNode.add(geneName);
            }
            associationNode.set("ensemblMappedGenes", mappedGeneNode);

            associationNode.set("traitName", JsonNodeFactory.instance.arrayNode().add(diseaseTrait));

            ArrayNode positionLinks = JsonNodeFactory.instance.arrayNode();
            for(JsonNode snp : associationObject.get("snps")){
                for(JsonNode location: snp.get("locations")){
                    String chName = location.get("chromosomeName").asText();
                    String chPosition  = location.get("chromosomePosition").asText();
                        String regionName = location.get("region").get("name").asText();
                        String mappedGene = chName + "|" + chPosition + "|" + regionName;
                        positionLinks.add(mappedGene);
                }
            }
            associationNode.set("efoLink", addEfoLinks(associationObject));
            associationNode.set("positionLinks", positionLinks);
            associationArrayNode.add(associationNode);

            JsonNode studyNode = studyMap.get(associationId);
            associationNode.put("author_s",
                    studyNode.get("publicationInfo").get("author").get("fullname").textValue());
            associationNode.put("publicationDate",
                    studyNode.get("publicationInfo").get("publicationDate").asText());
            associationNode.put("pubmedId",
                    studyNode.get("publicationInfo").get("pubmedId").textValue());
            associationNode.put("accessionId",
                    studyNode.get("accessionId").textValue());
            for(JsonNode ancestry: studyNode.get("ancestries")){

            }
        }
        return associationArrayNode;
    }

    private ArrayNode addEfoLinks(JsonNode associationObject){
        ArrayNode efoLinkNode = JsonNodeFactory.instance.arrayNode();
        for(JsonNode efoNode: associationObject.get("efoTraits")){
            String uri = efoNode.get("uri").textValue();
            String shortForm = efoNode.get("shortForm").textValue();
            String trait = efoNode.get("trait").textValue();
            efoLinkNode.add(trait + "|" + shortForm + "|" + uri + "|" + shortForm);
        }
        return efoLinkNode;
    }
}