import React, { useState, useEffect } from 'react';
import { getCookie} from '../../../utils/utilFunctions';
import { useNavigate } from 'react-router-dom';
import schema from '../../../schema/react-json-schema/Workflows/clusteringUsingRaceIDSchema.json';
import Form from 'react-jsonschema-form';
import styled from 'styled-components';
import { Container, Button } from "reactstrap";
import Switch from 'react-switch';
import Toggle from 'react-toggle';
import 'react-toggle/style.css';




export function ClusteringUsingRaceIDComponent(props) {
    
    let jwtToken = getCookie('jwtToken');
    const [formData, setFormData] = useState({});
    const navigate = useNavigate();


    const uiSchema = {
        "clustering": {
          "classNames": "category",
          "distanceMetric": {
            "classNames": "sub-category",
            "ui:title": "Distance Metric",
            "ui:widget": "select",
            "ui:options": {
              "placeholder": "Select a distance metric",
              "enumOptions": [
                { "label": "Euclidean", "value": "euclidean" },
                { "label": "Manhattan", "value": "manhattan" },
                { "label": "Cosine", "value": "cosine" }
              ]
            }
          },
          "clusteringMethod": {
            "ui:title": "Clustering Method",
            "ui:widget": "select",
            "classNames": "sub-category",
            "ui:options": {
              "placeholder": "Select a clustering method",
              "enumOptions": [
                { "label": "Hierarchical", "value": "hierarchical" },
                { "label": "K-Means", "value": "kmeans" },
                { "label": "DBSCAN", "value": "dbscan" }
              ]
            }
          },
          "default": {
            "ui:title": "Use Defaults?",
            "ui:widget": "select",
            "classNames": "category default-class",
            "ui:options": {
              "placeholder": "Use Defaults?",
              "classNames": "sub-category",
              "enumOptions": [
                { "label": "Yes", "value": "Yes" },
                { "label": "No", "value": "No" }
              ]
            }
          },
        },
        "outliers": {
          "classNames": "category",
          "minTranscripts": {
            "ui:title": "Minimum Transcripts",
            "classNames": "sub-category",
            "ui:widget": "updown",
            "ui:options": {
              "validate": false
            }
          },
          "minGenes": {
            "classNames": "sub-category",
            "ui:title": "Minimum Genes",
            "ui:widget": "updown",
            "ui:options": {
              "validate": false
            }
          },
          "plotFinalClusters": {
            "classNames": "sub-category",
            "ui:title": "Plot Final Clusters?",
            "ui:widget": "toggle"
          },
          "default": {
            "ui:title": "Use Defaults?",
            "ui:widget": "select",
            "classNames": "category default-class",
            "ui:options": {
              "placeholder": "Use Defaults?",
              "classNames": "sub-category",
              "enumOptions": [
                { "label": "Yes", "value": "Yes" },
                { "label": "No", "value": "No" }
              ]
            }
          },
        },
        "tSNEAndFR": {
          "classNames": "category",
          "perplexity": {
            "ui:title": "Perplexity",
            "classNames": "sub-category",
            "ui:widget": "updown",
            "ui:options": {
              "validate": false
            }
          },
          "KNN": {
            "classNames": "sub-category",
            "ui:title": "KNN",
            "ui:widget": "updown",
            "ui:options": {
              "validate": false
            }
          },
          "default": {
            "ui:title": "Use Defaults?",
            "ui:widget": "select",
            "classNames": "category default-class",
            "ui:options": {
              "placeholder": "Use Defaults?",
              "classNames": "sub-category",
              "enumOptions": [
                { "label": "Yes", "value": "Yes" },
                { "label": "No", "value": "No" }
              ]
            }
          },
        },
        "extraParameters": {
          "classNames": "category",
          "tableLimit": {
            "classNames": "sub-category",
            "ui:widget": "updown",
            "ui:options": {
              "validate": false
            }
          },
          "plotLimit": {
            "classNames": "sub-category",
            "ui:widget": "updown",
            "ui:options": {
              "validate": false
            }
          },
          "foldChange": {
            "classNames": "sub-category",
            "ui:widget": "updown",
            "ui:options": {
              "validate": false
            }
          },
          "pvalueCutoff": {
            "classNames": "sub-category",
            "ui:widget": "range",
          }
        }
      };

      const widgets = {
        toggle: (props) => (
          <Toggle
            checked={props.value}
            onChange={(e) => props.onChange(e.target.checked)}
          />
        ),
      };


  useEffect(() => {
    if(jwtToken) {

    } else {
    //   navigate("/routing");
    }
  }, [jwtToken, navigate]);

  return (
    <div className='workflow-container common-class-tools-and-workflows'>
      <div className="separator heading">
        <div className="stripe"></div> 
          <h2 className="h-sm font-weight-bold">
            Tool Parameters 
          </h2> 
        <div className="stripe"></div>
      </div>
        {/* <Container> */}
            <Form
                schema={schema}
                formData={formData}
                widgets={widgets}
                onChange={({ formData }) => setFormData(formData)}
                // onSubmit={handleSubmit}
                // SubmitButton={SubmitButton}
                uiSchema={uiSchema}
                // showErrorList={false}
                // noHtml5Validate={true}
            />
        {/* </Container> */}
    </div>
  )
};
